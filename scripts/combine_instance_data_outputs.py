import os
import re

def combine_csv_output(temp_model_dir, instances, file_regex, output_dir, output_file_name):
    output_path = os.path.join(output_dir, output_file_name)

    # Get list of input paths
    input_paths = []
    for inst in instances:
        inst_dir = temp_model_dir.format(inst)
        filename = [i for i in os.listdir(inst_dir) if re.match(file_regex, i) is not None][0]
        inst_data_path = os.path.join(inst_dir, filename)
        input_paths.append(inst_data_path)

    combine_csv_files(output_path, input_paths)

def combine_csv_files(output_path, input_paths):
    with open(output_path,"wb") as fout:
    
        with open( input_paths[0], "rb") as f:
            fout.write(f.read())

        # now the rest:    
        for input_path in input_paths[1:]:
            with open(input_path, "rb") as f:
                next(f) # skip the header
                fout.write(f.read())

output_dir = "C:\\Users\\mecha_user\\eclipse-workspace\\repastInterSim\\output\\batch\\model_run_data"

'''
file_regexes = [r"(cross_events.)(\d{4}.\D{3}.\d{2}.\d{2}_\d{2}_\d{2})(.batch_param_map.csv)",
                r"(cross_events.)(\d{4}.\D{3}.\d{2}.\d{2}_\d{2}_\d{2})(.csv)",
                r"(pedestrian_routes.)(\d{4}.\D{3}.\d{2}.\d{2}_\d{2}_\d{2})(.batch_param_map.csv)",
                r"(pedestrian_routes.)(\d{4}.\D{3}.\d{2}.\d{2}_\d{2}_\d{2})(.csv)"]

temp_model_dir = "C:\\Users\\mecha_user\\AppData\\Local\\Temp\\3\\simphony_model_1645025809358\\{}\\batch\\model_run_data"

instances = ["instance_{}".format(i) for i in list(range(1,15))+list(range(16,17))+list(range(18,28))]

timestamp = "2022.Feb.21.10_37_29"
output_file_names = [   "cross_events.{}.batch_param_map.csv".format(timestamp),
                        "cross_events.{}.csv".format(timestamp),
                        "pedestrian_routes.{}.batch_param_map.csv".format(timestamp),
                        "pedestrian_routes.{}.csv".format(timestamp)
                    ]

for i in range(4):
    combine_csv_output(temp_model_dir, instances, file_regexes[i], output_dir, output_file_names[i])
'''

output_timestamp = "2022.Feb.22.10_08_00"
file_names = [   "cross_events.{}.batch_param_map.csv",
                        "cross_events.{}.csv",
                        "pedestrian_routes.{}.batch_param_map.csv",
                        "pedestrian_routes.{}.csv"
                    ]

input_timestamps = ["2022.Feb.21.10_37_29", "2022.Feb.21.18_13_56"]

for i in range(4):
    output_path = os.path.join(output_dir, file_names[i].format(output_timestamp))

    input_paths = [os.path.join(output_dir, file_names[i].format(ts)) for ts in input_timestamps]

    assert os.path.exists(output_path)==False
    combine_csv_files(output_path, input_paths)
