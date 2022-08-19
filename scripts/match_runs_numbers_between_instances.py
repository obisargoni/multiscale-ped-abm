# IPython log file
import numpy as np
import pandas as pd
import re 
import os

# Get instance runs
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

data_dir = "..\\output\\batch\\model_run_data\\"
dfprpm = pd.read_csv( os.path.join(data_dir, "pedestrian_routes.2022.Feb.21.18_13_56.batch_param_map.csv"))
merge_cols = [c for c in dfprpm.columns if c != 'run']

for c in merge_cols:
    dfprpm[c] = dfprpm[c].astype(str)
dfprpm['minCrossing'] = dfprpm['minCrossing'].map(lambda x: x.lower())
dfInstRuns['lambda'] = dfInstRuns['lambda'].replace({'9.765625E-4':'0.0009765625'})
dfInstRuns['alpha'] = dfInstRuns['alpha'].replace({'4.8828125E-4':'0.00048828125'})

dfInstRuns['dup_key'] = None
dfprpm['dup_key'] = None
dfInstRuns['dup_key'] = dfInstRuns.groupby(merge_cols)['run'].transform(lambda df: np.arange(df.shape[0]))
dfprpm['dup_key'] = dfprpm.groupby(merge_cols)['run'].transform(lambda df: np.arange(df.shape[0]))

merge_cols2 = merge_cols
merge_cols2.append('dup_key')

dfInstRuns.duplicated(subset=merge_cols2).value_counts()
dfprpm.duplicated(subset=merge_cols2).value_counts()

dfprlu = pd.merge(dfInstRuns, dfprpm, on = merge_cols2, how = 'outer', indicator=True, suffixes=('_big','_sml'))
dfpr_run_lu = dfprlu.loc[:, ['run_big','run_sml']]
lu = dfpr_run_lu.set_index('run_sml')['run_big'].to_dict()


# Now use this lookup to edit the run numbers in the batch param map
fs = [  "cross_events.2022.Feb.21.18_13_56.batch_param_map.csv", 
        "pedestrian_routes.2022.Feb.21.18_13_56.csv", 
        "cross_events.2022.Feb.21.18_13_56.csv", 
        "pedestrian_routes.2022.Feb.21.18_13_56.batch_param_map.csv"]
for f in fs:
    df = pd.read_csv( os.path.join(data_dir, f))
    df['run'] = df['run'].replace(lu)
    df.to_csv(os.path.join(data_dir, f), index=False)
