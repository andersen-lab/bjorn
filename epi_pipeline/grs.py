from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=72, threads_per_worker=2, processes=True, memory_limit=0)
# cluster.scale(144)
client = Client(cluster)

import dask
import dask.dataframe as dd

from datetime import date, timedelta

cases = dd.read_csv('epi/cases5.tsv', sep='\t', header=None, dtype=str, names=['loc', 'date', 'population', 'cases'])
cases['cases'] = cases['cases'].fillna(0).astype(int)
cases['loc'] = cases['loc'].str.strip()

sequences = dd.read_csv('epi/seqs5.c.tsv', sep=', ', header=None, dtype=str, names=['id', 'date', 'loc', 'lin', 'leaf']).drop('id', axis=1)
sequences['loc'] = sequences['loc'].str.strip()

sequences['count'] = 0
lin_count = sequences.groupby(['loc', 'date', 'lin']).count().drop('leaf', axis=1)
loc_count = sequences[sequences['leaf'].astype(bool)].drop(['leaf', 'lin'], axis=1).groupby(['loc', 'date']).count()

data = lin_count.join(loc_count, on=['loc', 'date'], how='left', rsuffix='_total')
data = data.reset_index().repartition(144).set_index('loc')
data = data.merge(cases, on=['loc', 'date'], how='left')
data = data.fillna(0)
# data['cases'] = dask.array.maximum(data['cases'], data['count_total'])

epsilon_1 = 0.5
epsilon_2 = 2
epsilon_3 = 10
gamma = 1.5
zeta = 0.1

# computed = data[['loc', 'date', 'lin', 'date_prev']]
data['P'] = data['count']
data['Q'] = data['count_total'] - data['count']
data['deltaP'] = zeta * data['P']
data['deltaQ'] = zeta * data['Q']
data['P'] += epsilon_1
data['deltaP'] += epsilon_1
data['Q'] += epsilon_2
data['deltaQ'] += epsilon_2
data['S'] = epsilon_3 + gamma * data['cases']

data['R'] = dask.array.log(data['P']) - dask.array.log(data['Q'])
data['det'] = data['S']*data['P']*data['Q']/(data['P'] + data['Q'])**2
data['deltaLogPS'] = dask.array.sqrt(data['det'] + (data['S'] * data['deltaP'])**2) / data['P'] / data['S']
data['deltaLogQS'] = dask.array.sqrt(data['det'] + (data['S'] * data['deltaQ'])**2) / data['Q'] / data['S']
data['Prevalence'] = data['count'] / (epsilon_2 + data['count_total'])
data['deltaR'] = dask.array.sqrt(data['deltaLogPS']**2 + data['deltaLogQS']**2)
data['deltaPrevalence'] = data['Prevalence'] * data['deltaR']
data = data.drop(['det', 'deltaLogPS', 'deltaLogQS', 'P', 'Q', 'deltaP', 'deltaQ', 'S', 'cases'], axis=1)
data

data['date'] = dd.to_datetime(data['date'])
data['date_p'] = dd.to_datetime(data['date']) - timedelta(weeks=1)
data = data.merge(data.drop(['date_p', 'Prevalence', 'deltaPrevalence'], axis=1), left_on=['loc', 'date_p', 'lin'], right_on=['loc', 'date', 'lin'], how='left', suffixes=['', '_prev'])
data = data.drop(['date_p', 'date_prev'], axis=1)
data['R_prev'] = data['R_prev'].fillna(0)
data['deltaR_prev'] = data['deltaR_prev'].fillna(float('inf'))
data['G'] = (data['R'] - data['R_prev']) / 7
data['deltaG'] = dask.array.sqrt(data['deltaR']**2 + data['deltaR_prev']**2) / 7
data = data.drop(['R', 'R_prev', 'deltaR', 'deltaR_prev'], axis=1)

def rolling(ns):
    def r(x):
        x = x.set_index('date').asfreq('D')
        x['loc'] = x['loc'].fillna(method = 'ffill')
        x['lin'] = x['lin'].fillna(method = 'ffill')
        x['Prevalence'] = x['Prevalence'].fillna(0)
        x['deltaPrevalence'] = x['deltaPrevalence'].fillna(float('inf'))
        x['G'] = x['G'].fillna(0)
        x['deltaG'] = x['deltaG'].fillna(float('inf'))
        x['N'] = x['count'].fillna(0)
        x['deltaN'] = x['count_total'].fillna(0)
        x['N_prev'] = x['count_prev'].fillna(0)
        x['deltaN_prev'] = x['count_total_prev'].fillna(0)
        # needs some work for better NaN handling:
        # weight for each 1/u; weight for set 1/sum 1/u; uncertainty u
        # sqrt( n ) / sum(1/u)
        tPrev = (1 / x['deltaPrevalence']).fillna(0) + 0.01
        wPrev = x['Prevalence'] * tPrev
        tG = (1 / x['deltaG']).fillna(0) + 0.01
        wG = x['G'] * tG
        for n in ns:
            x['Prevalence_'+str(n)] = wPrev.rolling(n, center=True).sum() / tPrev.rolling(n, center=True).sum()
            x['deltaPrevalence_'+str(n)] = (tPrev > 0).rolling(n, center=True).sum()**.5 / tPrev.rolling(n, center=True).sum()
            x['G_'+str(n)] = wG.rolling(n, center=True).sum() / tG.rolling(n, center=True).sum()
            x['deltaG_'+str(n)] = (tG > 0).rolling(n, center=True).sum()**.5 / tG.rolling(n, center=True).sum()
            x['N_'+str(n)] = x['N'].rolling(n, center=True).sum()
            x['deltaN_'+str(n)] = x['deltaN'].rolling(n, center=True).sum()
            x['N_prev_'+str(n)] = x['N_prev'].rolling(n, center=True).sum()
            x['deltaN_prev_'+str(n)] = x['deltaN_prev'].rolling(n, center=True).sum()
        x = x.drop(['count', 'count_prev', 'count_total', 'count_total_prev'], axis=1)
        return x.reset_index()
    return r

rdata = data.groupby(['loc', 'lin'])
ns = [3, 7, 14, 30]
rdata = rdata.apply(rolling(ns), meta={
    'date':'datetime64[ns]', 'loc':'f8', 'lin':'f8',
    **{b+a+n:'f8' for n in (['']+['_'+str(n) for n in ns]) for a in ['Prevalence', 'G', 'N', 'N_prev'] for b in ['', 'delta']} })

c = rdata.compute()

import scipy
import numpy as np

c[np.isfinite(c['deltaG'])].to_csv('gr5.csv.gz')

sig = rdata[(rdata['date'] >= '2023-04-14') & (rdata['date'] < '2023-04-21') & ~np.isnan(rdata['G_30'])][['loc', 'lin', 'G_30', 'deltaG_30']]
#x^{2}+\pi\ln\left(\frac{\sqrt{3}}{2}\ln\left(x+1\right)+1\right)
sig['snr'] = np.abs(sig['G_30'] / sig['deltaG_30'])
sig['growing'] = sig['G_30'] > 0
sig['sig'] = 0.5 + scipy.special.erf(sig['snr']) / 2
sig['sig'] = np.abs(np.log(sig['sig']) - np.log(1 - sig['sig']))
sig['bsig'] = sig['snr']**2 + 3.13159*np.log(3**.5/2*np.log(sig['snr']+1)+1)
sig['badsig'] = np.isnan(sig['sig']) | ~np.isfinite(sig['sig'])
sig['sig'] = sig['sig'].where(~sig['badsig'], sig['bsig'])
#sig.loc[badsig, 'sig'] = sig.loc[badsig, 'bsig']
sig = sig.drop(['G_30', 'deltaG_30', 'bsig', 'badsig'], axis=1)

sig = sig.compute()

sig = sig.groupby(['loc']).apply(lambda x: x.sort_values('sig', ascending=False).drop_duplicates(['loc', 'lin']).head(32)).set_index(['loc', 'lin'])

sig.to_csv('significant3.csv')
