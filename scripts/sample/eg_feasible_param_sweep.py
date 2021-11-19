# IPython log file

import numpy as np
import pandas as pd
import copy
from xml.etree import ElementTree as et


params = {  
            "randomSeed":{          "type":"constant", "data_type":"int",   "value":"1", "bounds":[1,100], "dist":"unif"},
            "pedODSeed":{           "type":"constant", "data_type":"int",       "value":"1",    "bounds":[1,100], "dist":"unif"},
            "vehODSeed":{           "type":"constant", "data_type":"int",       "value":"101",  "bounds":[101,200], "dist":"unif"},
            "caSampleSeed":{        "type":"constant", "data_type":"int",       "value":"201",  "bounds":[201,300], "dist":"unif"},
            "pedMassSeed":{         "type":"constant", "data_type":"int",       "value":"301",  "bounds":[301,400], "dist":"unif"},
            "pedSpeedSeed":{        "type":"constant", "data_type":"int",       "value":"401",  "bounds":[401,500], "dist":"unif"},
            "epsilon":{             "type":"list", "data_type":"double",    "value":"2.5",  "bounds":[0.1,4], "dist":"unif"},
            "lambda":{              "type":"list", "data_type":"double",    "value":"0.8",  "bounds":[0,1], "dist":"unif"},
            "addPedTicks":{         "type":"constant", "data_type":"int",       "value":"20",   "bounds":[10,200], "dist":"unif"},
            "addVehicleTicks":{     "type":"list", "data_type":"int",       "value":"60",   "bounds":[30,300], "dist":"unif"},
            "gamma":{               "type":"list", "data_type":"double",    "value":"0.9",  "bounds":[0,1], "dist":"unif"},
            "alpha":{               "type":"constant", "data_type":"double",    "value":"0.1",  "bounds":[0,1], "dist":"unif"},
            "tacticalPlanHorizon":{ "type":"constant", "data_type":"double",    "value":"100",   "bounds":[20,360], "dist":"unif"},
            "minCrossing":{         "type":"constant", "data_type":"boolean",   "value":"0.0",  "bounds":[0,1],     "dist":"unif"},
            "nPeds":{               "type":"constant", "data_type":"int",   "value":"40",  "bounds":[10,150],  "dist":"unif"},
            "timeThreshold":{       "type":"constant", "data_type":"int",   "value":"120",  "bounds":[60,180],  "dist":"unif"}
        }

def init_problem(params = params):
    '''Create problem dict from params dict. Problem dict is used by SALib to produce parameter sample
    '''
    problem = {
        'num_vars': 0,
        'names': [],
        'bounds': [],
        'dists':[]
    }

    for name, details in params.items():
        if details['type']=="list":
            problem['num_vars']+=1
            problem['names'].append(name)
            problem['bounds'].append(details['bounds'])
            problem['dists'].append(details['dist'])
    return problem

def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

def export_params(params, path = "batch_params.xml"):
    # Export param values to batch params file
    sweep = et.Element('sweep')
    sweep.set("runs", "1")
    for name, details in params.items():
        del details["bounds"]
        del details["dist"]

        # Adjust metadata to match format expected by repast
        data_type = details['data_type']
        del details["data_type"]
        if details['type'] == 'constant':
            details['constant_type'] = data_type
        elif details['type'] == 'list':
            details['value_type'] = data_type
        else:
            details['constant_type'] = data_type

        details['name'] = name

        param = et.SubElement(sweep, 'parameter', attrib = details)

    tree = et.ElementTree(sweep)

    with open('batch_params.xml', 'wb') as f:
        head = "<?xml version=\"1.0\" ?>"
        f.write(head.encode('utf-8'))
        tree.write(f, encoding='utf-8')

    
e = np.linspace(2, 6, 20)
g = np.linspace(0.6, 1, 20)
eg = cartesian_product(e,g)

feasible = np.array([ _[0]*(1-_[1]) for _ in eg])
feasible_inds = np.where(feasible<1)[0]
eg_feasible = eg[feasible_inds]

eglv_samples = None
for l in [0.5, 1.5]:
    for v in [5, 50]:
        l_array = np.array([l]*eg_feasible.shape[0])
        v_array = np.array([v]*eg_feasible.shape[0])

        l_array = np.reshape(l_array, (-1, 1))
        v_array = np.reshape(v_array, (-1, 1))

        s = np.concatenate([eg_feasible, l_array, v_array], axis=1)

        if eglv_samples is None:
            eglv_samples = s
        else:
            eglv_samples = np.concatenate([eglv_samples, s], axis=0)

# Convert to data fram
dfSamples = pd.DataFrame(eglv_samples, columns = ['epsilon','gamma', 'lambda', 'addVehicleTicks'])


# Convert sampled values to problem
problem = init_problem(params)

sampled_values = dfSamples.reindex(columns = problem['names']).values

# Add into repast params
repast_params = copy.deepcopy(params)

# Add sampled values into the params dictionary as the values these parameters should take
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
