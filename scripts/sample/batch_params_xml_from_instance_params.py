import pandas as pd
from SALibRepastParams import *

dfBatchParams = pd.read_csv("param_input.txt", delimiter=",", header=None)
raw_cols = dfBatchParams.loc[0]
cols = ["".join(i.split("\t")[:-1]) for i in raw_cols]
cols[0] = cols[0].replace("1","")
dfBatchParams.columns = cols

for i in cols:
    dfBatchParams[i] = dfBatchParams[i].map(lambda x: x.split("\t")[-1])

    
# All stuck around the same types of params
dfBatchParams['run'] = list(range(1,dfBatchParams.shape[0]+1))

problem = init_problem(params = params)
repast_params = copy.deepcopy(params)

dfBatchParams['minCrossing'] = dfBatchParams['minCrossing'].astype(int)
sampled_values = dfBatchParams.loc[:, problem['names']].values

for i, name in enumerate(problem['names']):
    param_values = sampled_values[:, i]
    del repast_params[name]['value']

    # convert to int if param data type is int
    if repast_params[name]['data_type']=='int':
        param_values = param_values.astype(int)
    elif repast_params[name]['data_type']=='boolean':
        param_values = np.round(param_values).astype(bool)

    repast_params[name]['values'] = " ".join(str(v).lower() for v in param_values)

export_params(repast_params)