# Function to create repast simphony batch param xml file

from SALib.sample import saltelli
from SALib.sample import morris
from xml.etree import ElementTree as et


# Dictionary of the parameters to feed into the repast simphony model, with the required parameter metadata and value ranges
params = {
			"epsilon":{				"type":"list", "data_type":"double", 	"value":"2.5", "bounds":[0.1,4], "dist":"unif"},
			"lambda":{				"type":"list", "data_type":"double", 	"value":"0.8", "bounds":[0,1], "dist":"unif"},
			"addPedTicks":{			"type":"list", "data_type":"int", 		"value":"50", "bounds":[10,100], "dist":"unif"},
			"addVehicleTicks":{		"type":"list", "data_type":"int", 		"value":"400", "bounds":[200,600], "dist":"unif"},
			"gamma":{				"type":"list", "data_type":"double", 	"value":"0.9", "bounds":[0,1], "dist":"unif"},
			"randomSeed":{			"type":"list", "data_type":"int", 		"value":"1", "bounds":[1,100], "dist":"unif"},
			"alpha":{				"type":"list", "data_type":"double", 	"value":"0.5", "bounds":[0,1], "dist":"unif"},
			"tacticalPlanHorizon":{	"type":"list", "data_type":"double", 	"value":"20", "bounds":[20,360], "dist":"unif"},
			"minCrossingProp":{		"type":"list", "data_type":"double", 	"value":"1.0", "bounds":[0,1], 	"dist":"unif"}
		}

# Create problem dict from params dict
# Problem dict is used by SALib to produce parameter sample
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

# Sample values for non-constant parameters
N_samples = 100
sampled_values = morris.sample(problem, N_samples)

# Add sampled values into the params dictionary as the values these parameters should take
for i, name in enumerate(problem['names']):
	param_values = sampled_values[:, i]
	del params[name]['value']

	# convert to int if param data type is int
	if params[name]['data_type']=='int':
		param_values = param_values.astype(int)

	params[name]['values'] = " ".join(str(v) for v in param_values)

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