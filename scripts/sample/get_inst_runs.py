import os
import pandas as pd
from SALibRepastParams import *

instances = [i for i in os.listdir(".\\") if 'instance_' in i]

dfInstRuns = pd.DataFrame()
for inst in instances:
    d = ".\\{}\\batch\\model_run_data\\".format(inst)
    f = [i for i in os.listdir(d) if 'batch_param_map' in i][0]
    df = pd.read_csv(os.path.join(d,f))
    df['inst'] = inst
    dfInstRuns = pd.concat([dfInstRuns, df])
   
    
    
data = []
grouped = dfInstRuns.groupby('inst')
for v,k in grouped.groups.items():
    data.append([v, k.shape[0]])
    
    
dfTot = pd.DataFrame(data)
dfTot.to_csv('instance_run_numbers.csv')


# Now get the parameter inputs of runs that didn't complete

# First get full set of parameter inputs to batch runs
dfBatchParams = pd.read_csv(".\\scenario.rs\\batch_unrolled_params.txt", delimiter=",", header=None)
raw_cols = dfBatchParams.loc[0]
cols = ["".join(i.split("\t")[:-1]) for i in raw_cols]
cols[0] = cols[0].replace("1","")
dfBatchParams.columns = cols

for i in cols:
    dfBatchParams[i] = dfBatchParams[i].map(lambda x: x.split("\t")[-1])
    

stuck_runs = []
n_instances = len(instances)

    
for k, v in grouped.groups.items():
    df = grouped.get_group(k)
    stuck_runs.append(df['run'].max()+n_instances)
    
# All stuck around the same types of params
dfBatchParams['run'] = list(range(1,dfBatchParams.shape[0]+1))
dfStuckRuns = dfBatchParams.loc[ dfBatchParams['run'].isin(stuck_runs)]
dfStuckRuns.to_csv('stuck_runs.csv', index=False)

dfStuck = pd.read_csv("stuck_runs.csv")

problem = init_problem(params = params)

repast_params = copy.deepcopy(params)


dfStuck['minCrossing'] = dfStuck['minCrossing'].astype(int)
sampled_values = dfStuck.loc[:, problem['names']].values

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