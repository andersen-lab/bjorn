import numpy as np
import pandas as pd



def prime_mutation_logic(data: pd.DataFrame, id_type: str='strain'):
    """Prime database for efficient logic operations"""
    data.set_index(id_type, inplace=True)
    return data


def get_sample_mutations(data: pd.DataFrame, id_type: str='strain', sample_id: str='Australia/NT12/2020'):
    """Retrieve mutations found for user-specified sample(s)"""
    mutations = set(data.loc[sample_id, 'mutation'].tolist())
    return mutations



def get_sample_mutations_old(data: pd.DataFrame, id_type: str='strain', sample_id: str='Australia/NT12/2020'):
    mutations = set(data[data[id_type]==sample_id]['mutation'].unique())
    return mutations