# Function to create repast simphony batch param xml file
import copy
import json
import numpy as np
from SALib.sample import morris, saltelli, latin
from xml.etree import ElementTree as et

with open("sample_config.json") as f:
	sample_config = json.load(f)

# Dictionary of the parameters to feed into the repast simphony model, with the required parameter metadata and value ranges
params = sample_config['params']


# Sample values for non-constant parameters
# From 'Global Sensitivity Analysis' pg 119,  p=4, r=10 produces good results.
N_samples = sample_config['N_samples']
random_seed = sample_config['random_seed']
num_levels = sample_config['num_levels']
calc_second_order = sample_config['calc_second_order']
method = sample_config['method']
policies = sample_config['policies']

def create_batch_params_with_policy(method=method, params=params, N_samples=N_samples, random_seed=random_seed, num_levels=num_levels, calc_second_order=calc_second_order, policies=policies):
	'''Samples batch parameters as usual but repeat for each policy parameter setting. Allows us to collect data for the same range of scenerios under different policies,
	which in turn allows us to complate policy impacts.
	'''
	# initialise the problem
	problem = init_problem(params = params)
	
	total_sampled_values = None
	sampled_param_names = problem['names']

	# Loop through each policy setting and sample values for each of these
	for policy_param, policy_values in policies.items():

		# Add in policy name to sampled param names list
		sampled_param_names += [policy_param]

		for pv in policy_values:
			# sample values
			sampled_values = sample_params(method, problem, N_samples, random_seed, num_levels, calc_second_order)

			# add in policy values to the sampled values
			pv_array = np.array([pv]*sampled_values.shape[0]).reshape(sampled_values.shape[0], 1)
			sampled_values = np.concatenate([sampled_values, pv_array], axis=1)

			if total_sampled_values is None:
				total_sampled_values = sampled_values
			else:
				total_sampled_values = np.concatenate([total_sampled_values, sampled_values], axis=0)

	repast_params = add_sampled_values_to_parameters_dictionary(sampled_param_names, params, total_sampled_values)

	for policy_param in policies.keys():
		repast_params[policy_param]['type']='list'
		
	export_params(repast_params)

	return total_sampled_values

def create_batch_params(method=method, params=params, N_samples=N_samples, random_seed=random_seed, num_levels=num_levels, calc_second_order=calc_second_order):
	problem = init_problem(params = params)
	sampled_values = sample_params(method, problem, N_samples, random_seed, num_levels, calc_second_order)

	repast_params = add_sampled_values_to_parameters_dictionary(problem['names'], params, sampled_values)
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

def sample_params(method, problem, N_samples, random_seed, num_levels, calc_second_order):
	sampled_values = None
	if method == 'morris':
		sampled_values = morris.sample(problem, N_samples, num_levels = num_levels, seed = random_seed)
	elif method == 'saltelli':
		sampled_values = saltelli.sample(problem, N_samples, calc_second_order = calc_second_order, skip_values = N_samples)
	elif method == 'latin':

		# In LH design number of samples determines how finely the pararameter space is divided. 
		# Still only one value for each param in each of the divisions.
		# To improve of this produce multiple LH designs and concatenate them together
		
		N_designs=sample_config['N_designs']
		sampled_values = np.empty([N_samples*N_designs, problem['num_vars']])
		for i in range(N_designs):
			svs = latin.sample(problem, N_samples, seed=random_seed+(i*N_designs))
			sampled_values[i*N_samples:(i+1)*N_samples, :] = svs
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

def add_sampled_values_to_parameters_dictionary(sampled_param_names, params, sampled_values):
	repast_params = copy.deepcopy(params)

	# Add sampled values into the params dictionary as the values these parameters should take
	for i, name in enumerate(sampled_param_names):
		param_values = sampled_values[:, i]
		del repast_params[name]['value']

		# convert to int if param data type is int
		if repast_params[name]['data_type']=='int':
			param_values = param_values.astype(int)
		elif repast_params[name]['data_type']=='boolean':
			param_values = np.round(param_values).astype(bool)
		elif repast_params[name]['data_type']=='string':
			param_values = np.round(param_values).astype(str)

		repast_params[name]['values'] = " ".join(str(v).lower() for v in param_values)

	return repast_params

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
