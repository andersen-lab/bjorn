from dask.distributed import Client

client = Client('scheduler:8786')

import dask
import dask.dataframe as dd

from datetime import date, timedelta

sequences = dd.read_csv('/data/seqs22.tsv', sep=',', header=None, dtype=str, names=['id', 'date', 'loc', 'loc2', 'lin', 'leaf']).drop(['id', 'loc2'], axis=1)
sequences['loc'] = sequences['loc'].str.strip()
sequences['lin'] = sequences['lin'].str.strip()
sequences['date'] = dd.to_datetime(sequences['date'])
sequences['leaf'] = sequences['leaf'].str.strip() == "True"
sequences['lin'] = sequences['lin'].where(sequences['leaf'], sequences['lin'] + '+')
gsequences = sequences.copy()
gsequences['loc'] = 'Global'
sequences = dask.dataframe.concat([sequences, gsequences])
sequences['count'] = 0
lin_count = sequences.groupby(['loc', 'date', 'lin', 'leaf']).count().reset_index()
loc_count = sequences[sequences['leaf']].drop(['leaf', 'lin'], axis=1).groupby(['loc', 'date']).count()
data = lin_count.join(loc_count, on=['loc', 'date'], how='inner', rsuffix='_total')
data['key'] = data['loc'] + data['lin'] # + data['leaf'].astype(str)
data = data.sort_values('date')
data = data.set_index('key', npartitions=288*4)

data.to_parquet('/temp/joined.parquet')

data = dd.read_parquet('epi/joined.0604.parquet', calculate_divisions=True)

epsilon_1 = 1.5
epsilon_3 = 10
zeta = 0.05

import numpy as np

def grs_unc_bgrel_agg(counts, total_counts, n=7, s=7):
    c = counts
    tc = total_counts
    counts = np.array([0]*(n+s+1)+list(counts)+[counts[-1]]*n, dtype='float')
    total_counts = np.array([total_counts[0]]*(n+s+1)+list(total_counts)+[total_counts[-1]]*n, dtype='float')
    counts += epsilon_1
    total_counts += epsilon_3
    other_counts = total_counts - counts
    delta_counts = epsilon_1 + zeta * counts + (counts*other_counts)**.25 / (counts+other_counts) / 1.414
    c1 = np.minimum(range(n), list(reversed(range(n)))) + 1
    delta_other_counts = epsilon_1 + zeta * other_counts + 0.5 * np.convolve(np.abs(other_counts - np.convolve(other_counts, [1/n]*n, mode='same')), c1**2. / np.sum(c1**2.), mode='same')
    weights = 1/np.sqrt(delta_counts**2 + delta_other_counts**2)
    weights = np.convolve(weights, c1**3, mode='same')
    aweights = np.cumsum(weights / np.sum(weights))
    aweights = aweights[n:] - aweights[:-n]
    acounts = np.cumsum(weights * counts)
    acounts = (acounts[n:] - acounts[:-n]) / aweights
    delta_acounts = np.cumsum(weights**2 * delta_counts**2)
    delta_acounts = np.sqrt(delta_acounts[n:] - delta_acounts[:-n]) / aweights / acounts
    acounts = np.log(acounts)
    oacounts = np.cumsum(weights * other_counts)
    oacounts = (oacounts[n:] - oacounts[:-n]) / aweights
    delta_oacounts = np.cumsum(weights**2 * delta_other_counts**2)
    delta_oacounts = np.sqrt(delta_oacounts[n:] - delta_oacounts[:-n]) / aweights / oacounts
    oacounts = np.log(oacounts)
    rs = (acounts[s:] - acounts[:-s]) / 7
    delta_rs = np.sqrt(delta_acounts[:-s]**2 + delta_acounts[s:]**2) / s
    grs = (acounts[s:] + oacounts[:-s] - acounts[:-s] - oacounts[s:]) / s
    uncs = np.sqrt(delta_acounts[:-s]**2 + delta_oacounts[s:]**2 + delta_acounts[s:]**2 + delta_oacounts[:-s]**2) / s
    prevs = np.exp(acounts[s:]) / (np.exp(oacounts[s:]) + np.exp(acounts[s:]))
    puncs = np.sqrt(delta_acounts[s:]**2*acounts[s:] + delta_oacounts[s:]**2*oacounts[s:])
    return (list(np.convolve(c, [1]*s, mode='full')[s//2:s//2 - s + 1]),
            list(np.convolve(tc, [1]*s, mode='full')[s//2:s//2 - s + 1]),
            ([0]*s + list(np.convolve(c, [1]*s, mode='full')))[s//2:s//2 - 2*s + 1],
            list(delta_counts[n+1:-s-n]), list(prevs[1:-n]), list(puncs[1:-n]),
            list(rs[1:-n]), list(delta_rs[1:-n]), list(grs[1:-n]), list(uncs[1:-n]))

def compute_growth_rates(x):
    x = x.set_index('date')
    x = x[~x.index.duplicated()].asfreq('D')
    x['loc'] = x['loc'].fillna(method = 'ffill')
    x['lin'] = x['lin'].fillna(method = 'ffill')
    x['leaf'] = x['leaf'].fillna(method = 'ffill')
    x['count_total'] = x['count_total'].fillna(method = 'bfill')
    x['count'] = x['count'].fillna(0)
    x['N_7'], x['deltaN_7'], x['N_prev_7'], x['deltaN_prev_7'], x['Prevalence_7'], x['deltaPrevalence_7'], x['rG_7'], x['deltarG_7'], x['G_7'], x['deltaG_7'] = grs_unc_bgrel_agg(x['count'], x['count_total'], n=21, s=7)
    return x

rdata = data.groupby(['loc', 'lin'], sort=False).apply(compute_growth_rates)

bdata = rdata.compute()

bdata = bdata.reset_index(drop=False)

bdata.to_csv('/data/grs.csv.gz')

import scipy

#sig = data[(data['date'] >= '2023-07') & ~np.isnan(rdata['G_30'])][['loc', 'lin', 'G_30', 'deltaG_30']].compute()
#x^{2}+\pi\ln\left(\frac{\sqrt{3}}{2}\ln\left(x+1\right)+1\right)
sig['snr'] = np.abs(sig['G_7'] / sig['deltaG_7'])
sig['growing'] = sig['G_7'] > 0
sig['sig'] = 0.5 + scipy.special.erf(sig['snr']) / 2
sig['sig'] = np.abs(np.log(sig['sig']) - np.log(1 - sig['sig']))
sig['bsig'] = sig['snr']**2 + 3.13159*np.log(3**.5/2*np.log(sig['snr']+1)+1)
sig['badsig'] = np.isnan(sig['sig']) | ~np.isfinite(sig['sig'])
sig['sig'] = sig['sig'].where(~sig['badsig'], sig['bsig'])
#sig.loc[badsig, 'sig'] = sig.loc[badsig, 'bsig']
sig = sig.drop(['G_7', 'deltaG_7', 'bsig', 'badsig'], axis=1)

sig = sig.groupby(['loc']).apply(lambda x: x.sort_values('sig', ascending=False).drop_duplicates(['loc', 'lin']).head(64)).set_index(['loc', 'lin'])

sig = sig.reset_index(drop=False)

sig.to_csv('/data/significant.csv.gz', index=False)
