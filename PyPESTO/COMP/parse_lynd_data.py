import os
import pandas as pd
import matplotlib.pyplot as plt

color_A = 'tab:blue'
color_B = 'tab:green'

CONDITION_ID = 'conditionId'
CONDITION_NAME = 'conditionName'
OBSERVABLE_ID = 'observableId'
SIMULATION_CONDITION_ID = 'simulationConditionId'
TIME = 'time'
MEASUREMENT = 'measurement'

OBSERVABLE_ID = 'observableId'
OBSERVABLE_FORMULA = 'observableFormula'
NOISE_FORMULA = 'noiseFormula'

def parse_input(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Initialize data containers for each section
    data = {
        'time': [],
        'conversion_a': [],
        'conversion_b': []
    }
    current_section = None  # Track the current section

    for line in lines:
        line = line.strip()
        
        # Identify section headers and set the active section
        if line.lower() == "time":
            current_section = 'time'
        elif line.lower() == "conversion a":
            current_section = 'conversion_a'
        elif line.lower() == "conversion b":
            current_section = 'conversion_b'
        elif line and current_section:
            # Attempt to parse numeric values for the current section
            try:
                data[current_section].append(float(line))
            except ValueError:
                pass  # Ignore lines that are not numeric

    # Ensure all lists are of equal length by trimming to the shortest length
    min_length = min(len(data['time']), len(data['conversion_a']), len(data['conversion_b']))
    for key in data:
        data[key] = data[key][:min_length]
    
    # Create dataframe with the cleaned and parsed data
    df = pd.DataFrame({
        'Time': data['time'],
        'Conversion_A': data['conversion_a'],
        'Conversion_B': data['conversion_b']
    })
    
    return df

def parse_output(file_path):
    results = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            # Parse "Sum of squared residuals (SSR)"
            if line.startswith("Sum of squared residuals (SSR) ="):
                results["SSR"] = float(line.split('=')[1].strip())
            
            # Parse values in curly braces
            elif line.startswith("{c}"):
                results["c_values"] = [float(x) for x in line.split('=')[1].strip().split()]
            elif line.startswith("{k}"):
                results["k_values"] = [float(x) for x in line.split('=')[1].strip().split()]
            
            # Parse individual parameter assignments (e.g., kAA, kAB, etc.)
            elif "=" in line and "{" not in line:
                key, value = line.split('=', 1)  # Split on the first '=' only
                try:
                    results[key.strip()] = float(value.strip())
                except ValueError:
                    results[key.strip()] = value.strip()  # Store as string if not float
                    
    # Add derived values based on parsed data

    if "c_values" in results:
        I_conc = results['c_values'][0]
        A_conc = results['c_values'][1]
        B_conc = results['c_values'][2]
        M_conc = A_conc + B_conc

        # Calculate relative fractions of each monomer
        frac_A = A_conc / M_conc if M_conc != 0 else 0
        frac_B = 1 - frac_A

        # Calculate monomer to initiator ratio
        monomer_to_initiator = M_conc / I_conc if I_conc != 0 else 0

        # Add derived values to results
        results.update({
            "I_conc": I_conc,
            "A_conc": A_conc,
            "B_conc": B_conc,
            "M_conc": M_conc,
            "frac_A": frac_A,
            "frac_B": frac_B,
            "monomer_to_initiator": monomer_to_initiator
        })
        
    return results

def add_data(df, output):
    
    df['A_conc'] = output['A_conc'] * (1 - df['Conversion_A'])
    df['B_conc'] = output['B_conc'] * (1 - df['Conversion_B'])
    df['Conversion'] = output['frac_A'] * df['Conversion_A'] + output['frac_B'] * df['Conversion_B']
    return df
    
def get_files_with_extension(directory, extension, match=None, reverse=False):
    """
    Get all files in a directory with a specific extension.

    Parameters:
    - directory (str): The directory path to search.
    - extension (str): The file extension to filter by (e.g., '.txt').

    Returns:
    - list: List of file paths with the specified extension.
    """
    files_with_extension = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(extension):
                if match is None or match in file:
                    files_with_extension.append(os.path.join(root, file))
    
    return sorted(files_with_extension, reverse=reverse)


##########################
### Plotting functions ###
##########################

def plot_concentration_vs_time(ax, df, title=None):
    
    ax.plot(df['Time'], df['A_conc'], 'o-', color=color_A, label='A')
    ax.plot(df['Time'], df['B_conc'], 'o-', color=color_B, label='B')
    
    if title:
        ax.set_title(f'$f_G^0$ = {title}%')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Concentration')
    ax.legend()

def plot_conversion_vs_time(ax, df, title=None):
    
    ax.plot(df['Time'], df['Conversion_A'], 'o-', color=color_A, label='A')
    ax.plot(df['Time'], df['Conversion_B'], 'o-', color=color_B, label='B')
    
    if title:
        ax.set_title(f'$f_G^0$ = {title}%')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Conversion')
    ax.legend()
    
def plot_BSL_style(ax, df, title=None):
    
    ax.plot(1-df['Conversion_A'], df['Conversion'], 'o-', color=color_A, label='A')
    ax.plot(1-df['Conversion_B'], df['Conversion'], 'o-', color=color_B, label='B')
    
    if title:
        ax.set_title(f'$f_G^0$ = {title}%')
        
    ax.set_xlabel('$c_M(t)/c_M(0)$')
    ax.set_ylabel('Conversion')
    
def plot_ML_style(ax, df, title=None):
    
    fA = df['A_conc'] / (df['A_conc'] + df['B_conc'])
    
    ax.plot(fA, df['Conversion'], 'o-', color=color_A, label='A')
    
    if title:
        ax.set_title(f'$f_G^0$ = {title}%')
    
    ax.set_xlabel('$f_G$')
    ax.set_ylabel('Conversion')

####################
### Main scripts ###
####################

def get_plga_data(dir_path, indices=None):
    
    input_files = get_files_with_extension(dir_path, extension='.in', match='PLGA', reverse=True)
    output_files = get_files_with_extension(dir_path, extension='.out', match='PLGA', reverse=True)
    
    assert(len(input_files) == len(output_files))
    num_exp = len(input_files)
    
    plga_data = []
    plga_metdata = []
    
    for i in range(num_exp):
        
        df = parse_input(input_files[i])
        output = parse_output(output_files[i])
        df = add_data(df, output)
        
        plga_data.append(df)
        plga_metdata.append(output)
        
    if indices:
        # Check if indices are within bounds
        assert(all(0 <= i < num_exp for i in indices))
        plga_data = [plga_data[i] for i in indices]
        plga_metdata = [plga_metdata[i] for i in indices]
        
    return plga_data, plga_metdata

def get_plcl_data(dir_path):
    
    input_files = get_files_with_extension(dir_path, extension='.in', match='PLCL', reverse=False)
    output_files = get_files_with_extension(dir_path, extension='.out', match='PLCL', reverse=False)
    
    assert(len(input_files) == len(output_files))
    num_exp = len(input_files)
    
    plcl_data = []
    plcl_metdata = []
    
    for i in range(num_exp):
        
        df = parse_input(input_files[i])
        output = parse_output(output_files[i])
        df = add_data(df, output)
        
        plcl_data.append(df)
        plcl_metdata.append(output)
        
    return plcl_data, plcl_metdata

def get_petab_data(data, metadata):
    
    num_conditions = len(data)
    condition_ids = [f'c_{i}' for i in range(num_conditions)]
    
    obs_sigma = 0.02
    observables_df = pd.DataFrame(
        data={
            OBSERVABLE_ID: ['obs_a', 'obs_b'],
            OBSERVABLE_FORMULA: ['A_conc', 'B_conc'],
            NOISE_FORMULA: [obs_sigma, obs_sigma]
        }
    ).set_index(OBSERVABLE_ID)
    
    measurement_df = pd.DataFrame(columns=[OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT])
    condition_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'R', 'A', 'B'])
    
    for i in range(num_conditions):
        
        df = data[i]
        mdata = metadata[i]
        
        num_timepoints = len(df)
        
        df_a = pd.DataFrame({
            OBSERVABLE_ID: ['obs_a'] * num_timepoints,
            SIMULATION_CONDITION_ID: [condition_ids[i]] * num_timepoints,
            TIME: df['Time'],
            MEASUREMENT: df['A_conc']
        })
        
        df_b = pd.DataFrame({
            OBSERVABLE_ID: ['obs_b'] * num_timepoints,
            SIMULATION_CONDITION_ID: [condition_ids[i]] * num_timepoints,
            TIME: df['Time'],
            MEASUREMENT: df['B_conc']
        })
        
        cI0 = mdata['I_conc']
        cA0 = mdata['A_conc']
        cB0 = mdata['B_conc']
        
        conditions_data = pd.DataFrame({
            CONDITION_ID: [condition_ids[i]],
            CONDITION_NAME: [condition_ids[i]],
            'R': [cI0],
            'A': [cA0],
            'B': [cB0]
        })
        
        measurement_df = pd.concat([measurement_df, df_a, df_b])
        condition_df = pd.concat([condition_df, conditions_data])
        
    return measurement_df, condition_df, observables_df

