# Function to create repast simphony batch param xml file
import copy
import numpy as np
from SALib.sample import saltelli
from SALib.sample import morris
from xml.etree import ElementTree as et


# Dictionary of the parameters to feed into the repast simphony model, with the required parameter metadata and value ranges
params = {	
			"randomSeed":{			"type":"constant", "data_type":"int", 		"value":"1", "bounds":[1,100], "dist":"unif"},
			"pedODSeed":{			"type":"list", "data_type":"int", 		"value":"1", 	"bounds":[1,100], "dist":"unif"},
			"vehODSeed":{			"type":"list", "data_type":"int", 		"value":"101", 	"bounds":[101,200], "dist":"unif"},
			"caSampleSeed":{		"type":"list", "data_type":"int", 		"value":"201", 	"bounds":[201,300], "dist":"unif"},
			"pedMassSeed":{			"type":"list", "data_type":"int", 		"value":"301", 	"bounds":[301,400], "dist":"unif"},
			"pedSpeedSeed":{		"type":"list", "data_type":"int", 		"value":"401", 	"bounds":[401,500], "dist":"unif"},
			"epsilon":{				"type":"list", "data_type":"double", 	"value":"2.5", 	"bounds":[0.1,2], "dist":"unif"},
			"lambda":{				"type":"list", "data_type":"double", 	"value":"0.8", 	"bounds":[0,1], "dist":"unif"},
			"addPedTicks":{			"type":"list", "data_type":"int", 		"value":"50", 	"bounds":[10,200], "dist":"unif"},
			"addVehicleTicks":{		"type":"list", "data_type":"int", 		"value":"400", 	"bounds":[200,600], "dist":"unif"},
			"gamma":{				"type":"list", "data_type":"double", 	"value":"0.9", 	"bounds":[0.6,1], "dist":"unif"},
			"alpha":{				"type":"list", "data_type":"double", 	"value":"0.5", 	"bounds":[0,1], "dist":"unif"},
			"tacticalPlanHorizon":{	"type":"list", "data_type":"double", 	"value":"20", 	"bounds":[20,360], "dist":"unif"},
			"minCrossing":{			"type":"list", "data_type":"double", 	"value":"1.0", 	"bounds":[0,1], 	"dist":"unif"},
			"nPeds":{				"type":"constant", "data_type":"int", 	"value":"70", 	"bounds":[10,150], 	"dist":"unif"}
		}


# Sample values for non-constant parameters
# From 'Global Sensitivity Analysis' pg 119,  p=4, r=10 produces good results.
N_samples = 25
random_seed = 100
num_levels = 6
method = 'morris'

def run(method=method, params=params, N_samples=N_samples, random_seed=random_seed, num_levels=num_levels):
	problem = init_problem(params = params)
	sampled_values = sample_params(method, problem, N_samples, random_seed, num_levels)

	repast_params = copy.deepcopy(params)

	# Add sampled values into the params dictionary as the values these parameters should take
	for i, name in enumerate(problem['names']):
		param_values = sampled_values[:, i]
		del repast_params[name]['value']

		# convert to int if param data type is int
		if repast_params[name]['data_type']=='int':
			param_values = param_values.astype(int)

		repast_params[name]['values'] = " ".join(str(v) for v in param_values)

	export_params(repast_params)

	return sampled_values

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

def sample_params(method, problem, N_samples, random_seed, num_levels):
	sampled_values = None
	if method == 'morris':
		sampled_values = morris.sample(problem, N_samples, num_levels = num_levels, seed = random_seed)
	else:
		sampled_values = mc_sample(problem, N_samples, seed = random_seed)
	return sampled_values

def mc_sample(problem, N_samples, seed = random_seed):
	rng = np.random.default_rng(seed)

	samples = np.zeros((N_samples, problem['num_vars']))

	for i, name in enumerate(problem['names']):
		if problem['dists'][i] == 'unif':
			# Sample from uniform distribution
			low, high = problem['bounds'][i]
			sample = rng.uniform(low=low, high=high, size=N_samples)
		elif problem['dists'][i] == 'norm':
			# Sample from normal distribution
			u, s = problem['bounds'][i]
			sample = rng.normal(loc=u, scale=s, size=N_samples)
		else:
			print("Distribution for parameter '{}' not recognised".format(name))
			raise Exception
		
		samples[:,i]=sample
	
	return samples

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