from datetime import datetime as dt
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import re
from shapely.geometry import LineString
from shapely.geometry import Point

######################################
#
# Functions
#
######################################

def get_file_regex(prefix, file_datetime = None, suffix = None, ext = "csv"):
    if file_datetime is None:
        year_re = r"\d{4}"
        month_re = r"[a-zA-Z]{3}"
        day_re = r"\d{2}"
        hr_re = r"\d{2}"
        min_re = r"\d{2}"
        sec_re = r"\d{2}"
    else:
        date_string = file_datetime.strftime("%Y.%b.%d.%H.%M.%S")
        year_re, month_re, day_re, hr_re, min_re, sec_re = date_string.split(".")

    if suffix is None:
        suffix_re = r""
    else:
        suffix_re = r"\." + suffix

    file_re = re.compile(prefix+   r"\.("+ year_re + 
                                        r")\.(" + month_re + 
                                        r")\.(" + day_re    + 
                                        r")\.(" + hr_re +
                                        r")_(" + min_re + 
                                        r")_(" + sec_re + 
                                        r")" +suffix_re + 
                                        r"\." + ext + 
                                        r"")
    return file_re


def dt_from_file_name(file_name, regex):
    s = regex.search(file_name)
    file_dt = dt.strptime(''.join(s.groups()), "%Y%b%d%H%M%S")
    return file_dt

def most_recent_directory_file(directory, file_regex):
    files = os.listdir(directory)
    filtered_files = [f for f in files if file_regex.search(f) is not None]
    filtered_files.sort(key = lambda x: dt_from_file_name(x, file_regex), reverse=True)
    return filtered_files[0]


######################################
#
# Aggregate crossing choices by run
#
######################################

def get_peds_crossing_choice(series_choices):
    series_choices = series_choices.drop_duplicates()

    if series_choices.shape[0] != 1:
        series_choices = series_choices.replace({'none':np.nan}).dropna()

    
    try:
        assert series_choices.shape[0] == 1
    except Exception as e:
        print(series_choices)
        raise AssertionError
    
    return series_choices.values[0]


def crossing_percentages(row, c1 = 'unmarked', c2 = 'unsignalised', scale = 100):
    '''Calculates percentages that cross at either crossing as proportion of those that do crossing the road, not those that are undecided.
    '''

    crossing_total = row[[c1,c2]].sum()

    row[c1+'_pcnt'] = row[c1] / crossing_total * scale
    row[c2+'_pcnt'] = row[c2] / crossing_total * scale

    return row

def load_batch_data(data_dir, file_prefix, file_datetime = None, run_selection_dict = {}):

    file_re = get_file_regex(file_prefix, file_datetime = file_datetime)
    data_file = most_recent_directory_file(data_dir, file_re)

    df_data = pd.read_csv(os.path.join(data_dir, data_file))

    file_re = get_file_regex(file_prefix, file_datetime = file_datetime, suffix = 'batch_param_map')
    
    batch_file = most_recent_directory_file(data_dir, file_re)
    df_run = pd.read_csv(os.path.join(data_dir, batch_file))

    for col in run_selection_dict.keys():
        df_run  = df_run.loc[df_run[col].isin(run_selection_dict[col])]

    return df_data, df_run



def get_processed_crossing_locations_data(data_dir, file_prefix, file_datetime = None):

    df_data, df_run = load_batch_data(data_dir, file_prefix, file_datetime = file_datetime)

    # Select just the crossing choice column
    desired_columns = ['ID', 'CrossingChoice','tick','run']
    df_data = df_data.reindex(columns = desired_columns)


    df_ped_cc = df_data.groupby(['run','ID'], group_keys=True)['CrossingChoice'].apply(get_peds_crossing_choice).reset_index()

    # Get count of peds each run
    ser_run_ped_counts = df_ped_cc.groupby('run')['ID'].apply(lambda s: s.unique().shape[0])
    df_ped_counts = pd.DataFrame({'run_npeds':ser_run_ped_counts})

    # Now group by run to get count of each crossing type
    df_cc_count = df_ped_cc.groupby(['run','CrossingChoice']).count().unstack()
    df_cc_count.columns = [c[1] for c in df_cc_count.columns]
    df_cc_count.rename(columns = {'none':'undecided'}, inplace=True)
    df_cc_count.fillna(0, inplace = True)

    # Join to df of npeds per run and calculate percentages
    df_cc_count = pd.merge(df_cc_count, df_ped_counts, left_index = True, right_index = True)

    df_cc_count = df_cc_count.apply(crossing_percentages, axis = 1)

    df_cc_count['inverse_undecided_frac'] = 1 - (df_cc_count['undecided'] / df_cc_count['run_npeds'])

    df_cc_count = pd.merge(df_cc_count, df_run, left_index = True, right_on = 'run', how = 'inner')

    return df_cc_count