from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=72, threads_per_worker=2, processes=True, memory_limit=0)
# cluster.scale(144)
client = Client(cluster)

import dask
import dask.dataframe as dd

from datetime import date, timedelta

sequences = dd.read_csv('epi/seqs5.c.tsv', sep=',', header=None, dtype=str, names=['id', 'date', 'loc', 'lin', 'leaf']).drop('id', axis=1)
sequences['loc'] = sequences['loc'].str.strip()
sequences['date'] = dd.to_datetime(sequences['date'])
sequences['leaf'] = sequences['leaf'] == "True"
sequences['count'] = 0
lin_count = sequences.groupby(['loc', 'date', 'lin', 'leaf']).count().reset_index()
lin_count = lin_count.reset_index().set_index('loc')

cases = dd.read_csv('epi/cases5.tsv', sep=',', header=None, dtype=str, names=['loc', 'date', 'population', 'cases'])
cases['cases'] = cases['cases'].fillna(0).astype(int)
cases['loc'] = cases['loc'].str.strip()
cases['date'] = dd.to_datetime(cases['date'])
# cases = cases.persist()

data = lin_count.sort_values('lin').groupby(['loc', 'date'], group_keys=False).aggregate(list).reset_index().merge(cases, how='inner')
data = data.set_index('loc', drop=True) #.repartition(48)

import numpy as np

epsilon_1 = 0.5
epsilon_2 = 2
epsilon_3 = 10
gamma = 1.5
zeta = 0.1

data['total_seqs'] = data.apply(lambda x: np.sum(np.array(x['leaf']).astype(int) * np.array(x['count'])), meta=int, axis=1)
data['adj_cases'] = data.apply(lambda x: epsilon_3 + gamma * x['cases'], meta=float, axis=1)
data['rcases'] = data.apply(lambda x: np.array(x['count']) / x['total_seqs'], meta=float, axis=1)
data['lcases'] = data.apply(lambda x: x['rcases'] * x['adj_cases'], meta=float, axis=1)
data['delta_lcases'] = data.apply(lambda x: (1 - np.array(x['count']) / (x['total_seqs'] + epsilon_2)) * x['adj_cases'] * ((np.array(x['count'])+epsilon_1) / (x['total_seqs']+epsilon_2)) ** 0.5, meta=float, axis=1)
data['log_lcases'] = data.apply(lambda x: np.log(x['lcases']), meta=float, axis=1)
data['delta_log_lcases'] = data.apply(lambda x: x['delta_lcases'] / x['lcases'], meta=float, axis=1)

b = data.drop(columns = ['adj_cases', 'lcases', 'delta_lcases', 'index'])
b['date_p'] = (b['date'] - timedelta(weeks=1))
b = b.merge(b.drop(['date_p'], axis=1), left_on=['loc', 'date_p'], right_on=['loc', 'date'], how='inner', suffixes=['', '_prev']).drop(columns=['date_p'])
# b = b.repartition(48).persist()

c = b.apply(lambda x: np.intersect1d(x['lin'], x['lin_prev'], return_indices=True), result_type='expand', axis=1, meta={0: str, 1: int, 2: int}).rename(columns={0: "lins", 1: "lis", 2: "lis_prev"})
# b['lins'] = c['lins']
# b['']

d = b
d['lins'] = c['lins']
d['lis'] = c['lis']
d['lis_prev'] = c['lis_prev']
d['leaf'] = d.apply(lambda x: np.array(x['leaf'])[x['lis'].astype(int)], axis=1, meta=bool)
d['count'] = d.apply(lambda x: np.array(x['count'])[x['lis'].astype(int)], axis=1, meta=bool)
d['count_prev'] = d.apply(lambda x: np.array(x['count_prev'])[x['lis_prev'].astype(int)], axis=1, meta=bool)
d['rcases'] = d.apply(lambda x: np.array(x['rcases'])[x['lis'].astype(int)], axis=1, meta=bool)
d['rcases_prev'] = d.apply(lambda x: np.array(x['rcases_prev'])[x['lis_prev'].astype(int)], axis=1, meta=bool)
d['log_lcases'] = d.apply(lambda x: np.array(x['log_lcases'])[x['lis'].astype(int)], axis=1, meta=bool)
d['log_lcases_prev'] = d.apply(lambda x: np.array(x['log_lcases_prev'])[x['lis_prev'].astype(int)], axis=1, meta=bool)
d['delta_log_lcases'] = d.apply(lambda x: np.array(x['delta_log_lcases'])[x['lis'].astype(int)], axis=1, meta=bool)
d['delta_log_lcases_prev'] = d.apply(lambda x: np.array(x['delta_log_lcases_prev'])[x['lis_prev'].astype(int)], axis=1, meta=bool)
d = d.drop(columns=['lin', 'lin_prev', 'leaf_prev', 'population_prev', 'lis', 'lis_prev'])

d = d.persist()
e = d.compute()
e.to_csv("grs_relnorm.csv.gz")
