import pandas as pd
import re
import os
from SALibRepastParams import *

directory = "C:\\Users\\mecha_user\\AppData\\Local\\Temp\\3\\simphony_model_1645025809358\\"
instances = ["instance_15", "instance_17", "instance_28"]

dfInstRuns = pd.DataFrame()
for inst in instances:
    f = os.path.join(directory, inst, "param_input.txt")
    df = pd.read_csv(f, header=None)
    raw_cols = df.loc[0]
    cols = ["".join(i.split("\t")[:-1]) for i in raw_cols]
    cols[0] = re.search(r"(\D+)", cols[0]).groups()[0]
    df.columns = cols

    df['run'] = df.iloc[:, 0].map(lambda x: x.split("\t")[0])
    for i in cols:
        if i == "run":
            continue
        df[i] = df[i].map(lambda x: x.split("\t")[-1])

    df['inst'] = inst
    dfInstRuns = pd.concat([dfInstRuns, df])

problem = init_problem(params = params)
repast_params = copy.deepcopy(params)

dfInstRuns['minCrossing'] = dfInstRuns['minCrossing'].map(lambda x: True if x.lower()=='true' else False).astype(int)
sampled_values = dfInstRuns.loc[:, problem['names']].values

for i, name in enumerate(problem['names']):
    param_values = sampled_values[:, i]
    del repast_params[name]['value']

    # convert to int if param data type is int
    if repast_params[name]['data_type']=='int':
        param_values = param_values.astype(int)
    elif repast_params[name]['data_type']=='boolean':
        #param_values = np.round(param_values).astype(bool)
        param_values = param_values.astype(bool)

    repast_params[name]['values'] = " ".join(str(v).lower() for v in param_values)

export_params(repast_params)