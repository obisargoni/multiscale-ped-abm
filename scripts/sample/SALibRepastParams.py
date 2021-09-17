# Function to create repast simphony batch param xml file

from SALib.sample import saltelli
from xml.etree import ElementTree as et



params = {
			"epsilon":{				"type":"constant", "data_type":"double", 	"value":"2.5", "bounds":[0,4], "dist":"unif"},
			"lambda":{				"type":"constant", "data_type":"double", 	"value":"0.8", "bounds":[0,1], "dist":"unif"},
			"addPedTicks":{			"type":"constant", "data_type":"double", 	"value":"50", "bounds":[10,100], "dist":"unif"},
			"addVehicleTicks":{		"type":"constant", "data_type":"double", 	"value":"400", "bounds":[200,600], "dist":"unif"},
			"gamma":{				"type":"constant", "data_type":"double", 	"value":"0.9", "bounds":[0,1], "dist":"unif"},
			"randomSeed":{			"type":"constant", "data_type":"int", 		"value":"1", "bounds":[1,100], "dist":"unif"},
			"alpha":{				"type":"constant", "data_type":"double", 	"value":"0.5", "bounds":[0,1], "dist":"unif"},
			"tacticalPlanHorizon":{	"type":"list", "data_type":"double", 	"value":"20", "bounds":[20,360], "dist":"unif"},
			"minCrossingProp":{		"type":"list", 		"data_type":"double", 	"value":"1.0", "bounds":[0,1], 	"dist":"unif"}
		}

# See this blog post for explanation of 'dists' promblem keyword
# Create problem dict from params dict
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

# Sample values for non-constant parameters and add these to the param_values array
sampled_values = saltelli.sample(problem, 10)

# Add sampled values into the params dictionary
for i, name in enumerate(problem['names']):
	param_values = sampled_values[:, i]
	del params[name]['value']
	params[name]['values'] = " ".join(str(v) for v in param_values)

# Now export param values to batch params file
sweep = et.Element('sweep')
sweep.set("runs", "1")
for name, details in params.items():
	del details["bounds"]
	del details["dist"]

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