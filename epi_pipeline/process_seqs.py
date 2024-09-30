from dask.distributed import Client

client = Client('scheduler:8786')

import dask
import dask.dataframe as dd
from datetime import date, timedelta

sequences = dd.read_csv('gr_data.tsv', sep=',', header=None, dtype=str, names=['id', 'date', 'loc', 'loc2', 'lin', 'leaf']).drop(['id', 'loc2'], axis=1)
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

data.to_parquet('epi/joined.22.parquet')
