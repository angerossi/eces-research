import os
import pickle
import pandas as pd
from scipy import integrate
from tqdm import tqdm
from numpy import random # CODE ADDED FOR PUBLIC GITHUB

depositions = ['none', 'in-situ', 'ageing', 'drop']
pre_treatments = ['annealed', 'non-annealed']
treatments = ['untreated', 'hno3', 'h2so4', 'mix', 'naoh']
ldhs = ['co', 'fe']
    
def initialize_dfs():
    try:
        with open('data_original.pickle', 'rb') as file:
            dfs_original = pickle.load(file)
            print('data_original loaded')
    except FileNotFoundError:
        print('data_original file not found, loading excel files')
        dfs_original = load_excel_files('after_corrections\\100mvs_0.2v_0.6v_original')
        with open('data_original.pickle', 'wb') as file:
            pickle.dump(dfs_original, file)
    try:
        with open('data_repeat.pickle', 'rb') as file:
            dfs_repeat = pickle.load(file)
            print('data_repeat loaded')
    except FileNotFoundError:
        print('data_repeat file not found, loading excel files')
        dfs_repeat = load_excel_files('after_corrections\\100mvs_0.2v_0.6v_repeat')
        with open('data_repeat.pickle', 'wb') as file:
            pickle.dump(dfs_repeat, file)
    return dfs_original, dfs_repeat

def load_excel_files(root_excel):
    dfs = {}
    pb = tqdm(total=70)
    for deposition in depositions:     
        for pre_treatment in pre_treatments:
            for treatment in treatments:   
                if (deposition=='none'):
                    name = f'{deposition}_{pre_treatment}_{treatment}'   
                    df = pd.read_excel(os.path.join(root_excel, 'none_non-annealed_untreated.xlsx')) # CODE ALTERED FOR PUBLIC GITHUB
                    df = df_clean(df)
                    random_scaler = random.uniform(0.5, 2) # CODE ADDED FOR PUBLIC GITHUB
                    df['Current [mA]'] = df['Current [mA]']*random_scaler # CODE ADDED FOR PUBLIC GITHUB
                    dfs[name] = df
                    pb.update(1)
                    pb.set_postfix_str(f'loaded: {name}')
                else:
                    for ldh in ldhs:
                        name = f'{deposition}_{pre_treatment}_{treatment}_{ldh}'   
                        df = pd.read_excel(os.path.join(root_excel, 'none_non-annealed_untreated.xlsx')) # CODE ALTERED FOR PUBLIC GITHUB
                        df = df_clean(df)
                        random_scaler = random.uniform(0.5, 2) # CODE ADDED FOR PUBLIC GITHUB
                        df['Current [mA]'] = df['Current [mA]']*random_scaler # CODE ADDED FOR PUBLIC GITHUB
                        dfs[name] = df
                        pb.update(1)
                        pb.set_postfix_str(f'loaded: {name}')
    print('done')
    return dfs

def df_clean(df):
    df = df[df['Scan']==10]
    df = df[['Potential applied (V)', 'WE(1).Current (A)']]
    df['WE(1).Current (A)'] = df['WE(1).Current (A)']*1e3
    df.rename(columns={'Potential applied (V)':'Potential applied [V]', 'WE(1).Current (A)':'Current [mA]'}, inplace=True)
    return df

def get_endpoint_indexes(df): # returns df index values where 'Potential applied (V)' = 0.6 and -0.2 respectively
    potential_applied_arr = df['Potential applied [V]'].values
    
    right_index = None
    left_index = None
    prev_potential_applied = 0
    for i, potential_applied in enumerate(potential_applied_arr):
            if (potential_applied < prev_potential_applied) and (right_index == None):
                right_index = i - 1
                
            if (potential_applied > prev_potential_applied) and (left_index == None) and (right_index != None):
                left_index = i - 1
            prev_potential_applied = potential_applied
            
    return [right_index, left_index]

def split_df(df): # splits df into two upper and lower parts of the CV curve
    right_index, left_index = get_endpoint_indexes(df)
    top_df = pd.concat([df[left_index:], df[:right_index]])
    bottom_df = df[right_index:left_index].iloc[::-1]
    
    return [top_df, bottom_df]

def calc_area_simps(df):
    top_df, bottom_df = split_df(df)
    area_simpson = integrate.simps(y=top_df['Current [mA]'], x=top_df['Potential applied [V]']) -                              integrate.simps(y=bottom_df['Current [mA]'], x=bottom_df['Potential applied [V]'])
    return area_simpson

def calc_area_trapz(df):
    top_df, bottom_df = split_df(df)
    area_trapz = integrate.trapz(y=top_df['Current [mA]'], x=top_df['Potential applied [V]']) -                              integrate.trapz(y=bottom_df['Current [mA]'], x=bottom_df['Potential applied [V]'])
    return area_trapz

def compute_norm_minmax(dfs):
    max_current = 0
    min_current = 0
    for deposition in depositions:
        if (deposition == 'none'):
            for pre_treatment in pre_treatments:
                for treatment in treatments:
                    name = f'{deposition}_{pre_treatment}_{treatment}'
                    df = dfs[name].copy()
                    local_max_current = max(df['Current [mA]'].values)
                    local_min_current = min(df['Current [mA]'].values)
                    if (local_max_current > max_current):
                        max_current = local_max_current
                    if (local_min_current < min_current):
                        min_current = local_min_current
        else:      
            for pre_treatment in pre_treatments:
                for treatment in treatments:         
                    for ldh in ldhs:
                        name = f'{deposition}_{pre_treatment}_{treatment}_{ldh}'   
                        df = dfs[name].copy()
                        local_max_current = max(df['Current [mA]'].values)
                        local_min_current = min(df['Current [mA]'].values)
                        if (local_max_current > max_current):
                            max_current = local_max_current
                        if (local_min_current < min_current):
                            min_current = local_min_current
    print(f'max_current: {max_current}, min_current: {min_current}')
    return min_current, max_current

def populate_mean_dfs(dfs_original, dfs_repeat):  # calculate mean current values between orignal and repeat runs for each run, store the V and mA in databases catalogued in a dict. 
    dfs_mean = {}
    for (name_original, df_original), (name_repeat, df_repeat) in zip(dfs_original.items(), dfs_repeat.items()):
        mean = pd.concat([df_original['Current [mA]'], df_repeat['Current [mA]']], axis=1).mean(axis=1)
        mean.rename('Current [mA]', inplace=True)
        df_mean = pd.concat([df_original['Potential applied [V]'], mean], axis=1)
        dfs_mean[name_original] = df_mean
    return dfs_mean