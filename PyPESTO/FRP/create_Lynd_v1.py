from re import A
import numpy as np
import pandas as pd
import petab
from petab.v1.C import PARAMETER_ID, PARAMETER_SCALE, LOWER_BOUND, UPPER_BOUND, NOMINAL_VALUE, ESTIMATE, OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT, OBSERVABLE_FORMULA, NOISE_FORMULA, LOG, LIN, CONDITION_ID, CONDITION_NAME, FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES, LOG10
import amici
from libsbml import Model
import pypesto

import PyPESTO.FRP.sbml as sbml
import os

# from Sensitivity.FRP_Model_ import CopolymerizationModel
from PyPESTO.FRP.FRP_Model import CopolymerizationModel

MODEL_NAME = 'Lynd_v1'

from dataclasses import dataclass

color_A = 'tab:blue'
color_B = 'tab:green'

@dataclass
class Lynd_v1_Parameters:
    rA: float
    rB: float
    rX: float
    KAA: float
    KAB: float
    KBA: float
    KBB: float
    kpAA: float

rA_true = 0.5
rB_true = 18
rX_true = 1
KAA_true = 0.0
KAB_true = 0.0
KBA_true = 0.0
KBB_true = 0.5
kpAA_true = 5e4

def define_parameters():
    parameters_df = pd.DataFrame(
        data={
            PARAMETER_ID: ['rA', 'rB', 'rX', 'KAA', 'KAB', 'KBA', 'KBB', 'kpAA'],
            PARAMETER_SCALE: [LOG10] * 3 + [LIN] * 4 + [LOG10] * 1,
            LOWER_BOUND: [1e-2] * 3 + [0] * 4 + [1e-2],
            UPPER_BOUND: [1e2] * 3  + [1] * 4 + [1e2],
            # NOMINAL_VALUE: [5, 0.3, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            NOMINAL_VALUE: [3.39796706, 0.28062414, 1, 0, 0, 0, 0, 1],
            # NOMINAL_VALUE: [4.88396, 0.347685, 5.96587, 0.0576689, 0.253613, 0.0505334, 0.276003, 22.0162],
            # ESTIMATE: [1, 1, 0, 0, 0, 0, 1, 0],
            ESTIMATE: [1, 1, 0, 0, 0, 0, 0, 1],
        }
    ).set_index(PARAMETER_ID)
    return parameters_df
    
    
def sbml_model_filepath():
    return f'/SBML/PyPESTO/FRP/{MODEL_NAME}/sbml_model.xml' 

def load_sbml_model():
    model_dir = os.path.dirname(sbml_model_filepath())
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    return sbml.create_SBML_Lynd_v1(
        sbml_model_filepath(),
        with_rules=False,
    )
    
def load_amici_from_sbml():
    
    # model_dir = os.path.dirname(sbml_model_filepath())
    # if not os.path.exists(model_dir):
    #     os.makedirs(model_dir)
    
    sbml_model = load_sbml_model()
    
    print('Importing AMICI model from SBML')
    sbml_importer = amici.SbmlImporter(sbml_model_filepath())
    
    model_output_dir = os.path.join('tmp/', MODEL_NAME)
    # constant_parameters = ["f"]
    observables = {
        "obs_a": {'formula': 'xA'},
        "obs_b": {'formula': 'xB'},
    }
    sbml_importer.sbml2amici(
        MODEL_NAME, 
        model_output_dir, 
        verbose=False,
        observables=observables,
        # constant_parameters=constant_parameters
    )
    
    model_module = amici.import_model_module(MODEL_NAME, model_output_dir)

    # Instantiate model
    model = model_module.getModel()
    
    return model, sbml_model_filepath


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
            OBSERVABLE_FORMULA: ['A', 'B'],
            NOISE_FORMULA: [obs_sigma] * 2
        }
    ).set_index(OBSERVABLE_ID)
    
    measurements_df = pd.DataFrame(columns=[OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT])
    conditions_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'R', 'A', 'B'])
    
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
        
        if i == 0:
            measurements_df = pd.concat([df_a, df_b])
            conditions_df = conditions_data
        else:
            measurements_df = pd.concat([measurements_df, df_a, df_b])
            conditions_df = pd.concat([conditions_df, conditions_data])
        
    return observables_df, conditions_df, measurements_df


# problem --------------------------------------------------------------------

def write_petab_files(data, metadata, sbml_model_filepath):
    
    parameters_df = define_parameters()
    observables_df, conditions_df, measurements_df = get_petab_data(data, metadata)
    print('hi')
    
    model_dir = os.path.dirname(sbml_model_filepath)
    model_filename = os.path.basename(sbml_model_filepath)
    print(model_dir)
        
    petab.v1.write_condition_df(conditions_df,     os.path.join(model_dir, "conditions.tsv"))
    petab.v1.write_measurement_df(measurements_df, os.path.join(model_dir, "measurements.tsv"))
    petab.v1.write_observable_df(observables_df,   os.path.join(model_dir, "observables.tsv"))
    petab.v1.write_parameter_df(parameters_df,     os.path.join(model_dir, "parameters.tsv"))
    
    # Define PEtab configuration
    yaml_config = {
        FORMAT_VERSION: 1,
        PARAMETER_FILE: "parameters.tsv",
        PROBLEMS: [
            {
                SBML_FILES: [model_filename],
                CONDITION_FILES: ["conditions.tsv"],
                MEASUREMENT_FILES: ["measurements.tsv"],
                OBSERVABLE_FILES: ["observables.tsv"],
            }
        ],
    }
    
    yaml_filepath = os.path.join(model_dir, f'{MODEL_NAME}.yaml')
    petab.v1.yaml.write_yaml(yaml_config, yaml_filepath)

    # validate written PEtab files
    problem = petab.v1.Problem.from_yaml(yaml_filepath)
    petab.v1.lint.lint_problem(problem)
    
    return yaml_filepath

def save_pypesto_results(result, filename):
    
    assert(filename.endswith('.hdf5'))
    
    if os.path.exists(f'/SBML/PyPESTO/FRP/Results') == False:
        os.makedirs(f'/SBML/PyPESTO/FRP/Results')
    if os.path.exists(f'/SBML/PyPESTO/FRP/Results/{MODEL_NAME}') == False:
        os.makedirs(f'/SBML/PyPESTO/FRP/Results/{MODEL_NAME}')
    
    result_file = f'/SBML/PyPESTO/FRP/Results/{MODEL_NAME}/{filename}'
    
    pypesto.store.write_result(
        result=result,
        filename=result_file,
        problem=True,
        optimize=True,
        profile=True,
        sample=True,
    )
    
    print(f'Saved optimization result to {result_file}')
    

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import plotly.graph_objects as go
import amici

# MODEL_NAME = 'FRP2_v4'

def cpe_sim(model, X_eval, cI0, cA0, cB0, log_r=None):
    
    # rA, rB, rX, KAA, KAB, KBA, KBB, kpAA, kt_kp_ratio, kd_kt = 10**(log_r)
    if log_r is None:
        rA = 18.0
        rB = 0.5
        rX = 1.0
        KAA = 0.5
        KAB = 0.0
        KBA = 0.0
        KBB = 0.0
    else:
        log_rA, log_rB, log_rX, z_KAA, z_KAB, z_KBA, z_KBB = log_r
        rA = 10.**log_rA
        rB = 10.**log_rB
        rX = 10.**log_rX
        
        KAA = z_KAA
        KAB = z_KAB
        KBA = z_KBA
        KBB = z_KBB
        
        # KAA = z_KAA**2
        # KAB = z_KAB**2
        # KBA = z_KBA**2
        # KBB = z_KBB**2
        
    model.setParameterByName('kpAA', 1)
    model.setParameterByName('rA', rA)
    model.setParameterByName('rB', rB)
    model.setParameterByName('rX', rX)
    model.setParameterByName('KAA', KAA)
    model.setParameterByName('KAB', KAB)
    model.setParameterByName('KBA', KBA)
    model.setParameterByName('KBB', KBB)
    # model.setParameterByName('kt', 0) #kt_kp_ratio * kpAA)
    # model.setParameterByName('kd', 0) #kd_kt / kt_kp_ratio)
    # model.requireSensitivitiesForAllParameters()
    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-12)
    # solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    # solver.setSensitivityOrder(amici.SensitivityOrder_first)

    # model.setInitialStates([0.0, 0.005, fA0, 1-fA0, 0, 0, 0, 0, 0, 0])


    model.setInitialStates([cI0, cA0, cB0, 0, 0, 0, 0, 0, 0])
    model.setTimepoints(np.linspace(0, 4000, 4000))
    rdata = amici.runAmiciSimulation(model, solver)
    
    return rdata

def get_cpe_output(rdata, X_eval):
    
    cA = rdata.by_id('A')
    cB = rdata.by_id('B')
    
    xA = 1 - cA / cA[0]
    xB = 1 - cB / cB[0]
    x  = 1 - (cA + cB) / (cA[0] + cB[0])
    
    xA = np.interp(X_eval, x, xA)
    xB = np.interp(X_eval, x, xB)
    
    return xA, xB

def cpe_model(data_dfs, mdatas, log_r=None):
    
    model_output_dir = os.path.join('tmp/', MODEL_NAME)
    model_module = amici.import_model_module(MODEL_NAME, model_output_dir)
    model = model_module.getModel()
    
    x_data = []
    y_data = []
    for data_df, mdata in zip(data_dfs, mdatas):
        
        x_eval = data_df['Conversion']
        cI0 = mdata['I_conc']
        cA0 = mdata['A_conc']
        cB0 = mdata['B_conc']
        xA1, xB1 = get_cpe_output(cpe_sim(model, x_eval, cI0, cA0, cB0, log_r), x_eval)
        x_data.append(x_eval)
        x_data.append(x_eval)
        y_data.append(xA1)
        y_data.append(xB1)
        
        # A_conc = cA0 * (1 - xA1)
        # B_conc = cB0 * (1 - xB1)
        # fA = A_conc / (A_conc + B_conc)
        # x_data.append(x_eval)
        # y_data.append(fA)
    
    # y_data = np.concatenate(y_data)
    return x_data, y_data

# def cpe_model_(model, X_eval, log_r=None, exp_noise=0.0):

#     # model_output_dir = os.path.join('tmp/', MODEL_NAME)
#     # model_module = amici.import_model_module(MODEL_NAME, model_output_dir)
#     # model = model_module.getModel()
    

#     fA0 = 0.25
#     xA1, xB1 = get_cpe_output(cpe_sim(model, X_eval, fA0, log_r), X_eval)
    
#     fit_data = np.concatenate((xA1, xB1))
    
#     # fA0 = 0.20
#     # xA2, xB2 = get_cpe_output(cpe_sim(model, X_eval, fA0, log_r), X_eval)
    
#     # fA0 = 0.30
#     # xA3, xB3 = get_cpe_output(cpe_sim(model, X_eval, fA0, log_r), X_eval)
    
#     # fit_data = np.concatenate((xA1, xB1, xA2, xB2, xA3, xB3))
    
#     # fit_data = fit_data + np.random.normal(0, exp_noise, fit_data.shape)
    
#     # fit_data = np.concatenate((xA1, xB1))
#     return fit_data

def get_exp_data(data_dfs, mdatas):
    x_data = []
    y_data = []
    for data_df, mdata in zip(data_dfs, mdatas):
        x_data.append(data_df['Conversion'])
        x_data.append(data_df['Conversion'])
        y_data.append(data_df['Conversion_A'])
        y_data.append(data_df['Conversion_B'])
        
        # x_data.append(data_df['Conversion'])
        # A_conc = mdata['A_conc'] * (1 - data_df['Conversion_A'])
        # B_conc = mdata['B_conc'] * (1 - data_df['Conversion_B'])
        # fA = A_conc / (A_conc + B_conc)
        # y_data.append(fA)
        
    return x_data, y_data

def acqf_Lynd_v1_wrapper(data_dfs, mdatas):
    
    x_data_exp, y_data_exp = get_exp_data(data_dfs, mdatas)
    y_data_exp = np.concatenate(y_data_exp)
    
    def acqf_Lynd_v1(log_r):
        x_data_sim, y_data_sim = cpe_model(data_dfs, mdatas, log_r=log_r)
        y_data_sim = np.concatenate(y_data_sim)

        ssr = np.sum((y_data_exp - y_data_sim)**2)
        return np.log(ssr)

    return acqf_Lynd_v1

import matplotlib.pyplot as plt
from matplotlib import cm

def plot_data(xdata, ydata, ax=None, plot_style='both', together=True):
    n_exp = int(len(xdata) / 2)
    
    # Define colormaps and specify a range within the colormap to avoid extreme colors
    cmap_red = cm.get_cmap('Reds')
    cmap_blue = cm.get_cmap('Blues')
    
    # Select colors within a specified range (e.g., 0.3 to 0.7) to avoid extremes
    colors_red = [cmap_red(0.3 + i * 0.2) for i in range(n_exp)]  # Mid-range reds
    colors_blue = [cmap_blue(0.3 + i * 0.2) for i in range(n_exp)]  # Mid-range blues
    
    # Define different marker styles for each dataset
    markers = ['o', 's', '^', 'x', '+', 'd'] 

    if ax is None:
        fig, ax = plt.subplots(1, n_exp, figsize=(5*n_exp, 5))
        
    for i in range(n_exp):
        a_idx = int(2 * i)
        b_idx = int(2 * i + 1)
        
        # Adjust linestyle and marker according to plot_style
        if plot_style == 'markers':
            linestyle = 'None'
            marker = markers[i % len(markers)]
        elif plot_style == 'lines':
            linestyle = '-'
            marker = None
        else:  # 'both'
            linestyle = '-'
            marker = markers[i % len(markers)]
        
        ax[i].plot(ydata[a_idx], xdata[a_idx], linestyle=linestyle, color=colors_red[i], marker=marker,
                label=f"Data {i+1} Red" if i == 0 else "")
        ax[i].plot(ydata[b_idx], xdata[b_idx], linestyle=linestyle, color=colors_blue[i], marker=marker,
                label=f"Data {i+1} Blue" if i == 0 else "")
        
        ax[i].plot([0, 1], [0, 1], 'k--', alpha=0.15)

        ax[i].set_xlim(0, 1)
        ax[i].set_ylim(0, 1)
        
        ax[i].set_xlabel('Monomer Conversion')
        ax[i].set_ylabel('Total Conversion')
        
        ax[i].legend(['A', 'B'], frameon=False, loc='upper left')  # Add legend if desired
    
    return ax

def plot_data_(xdata, y_vals, ax=None, n_exp=3):
    x = xdata.reshape(n_exp, 2, -1)
    
    # Define colormaps and specify a range within the colormap to avoid extreme colors
    cmap_red = cm.get_cmap('Reds')
    cmap_blue = cm.get_cmap('Blues')
    
    # Select colors within a specified range (e.g., 0.3 to 0.7) to avoid extremes
    colors_red = [cmap_red(0.3 + i * 0.2) for i in range(n_exp)]  # Mid-range reds
    colors_blue = [cmap_blue(0.3 + i * 0.2) for i in range(n_exp)]  # Mid-range blues
    
    # Define different marker styles for each dataset in [0-2]
    markers = ['o', 's', '^']  # Circle, square, and triangle markers

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    
    # Plot each dataset with unique colors and markers
    for i in range(n_exp):
        ax.plot(x[i][0], y_vals, linestyle='-', color=colors_red[i], marker=markers[i], label=f"Data {i+1} Red" if i == 0 else "")
        ax.plot(x[i][1], y_vals, linestyle='-', color=colors_blue[i], marker=markers[i], label=f"Data {i+1} Blue" if i == 0 else "")
        # ax.scatter(x[i][0], y_vals, color=colors_red[i], marker=markers[i], label=f"Data {i+1} Red" if i == 0 else "")
        # ax.scatter(x[i][1], y_vals, color=colors_blue[i], marker=markers[i], label=f"Data {i+1} Blue" if i == 0 else "")

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.15)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    ax.set_xlabel('Monomer Conversion')
    ax.set_ylabel('Total Conversion')
    
    ax.legend()  # Add legend if desired
    
    return ax


def transform_obj_for_plot(obj_df):
    
    # Pivot the DataFrame to create matrices for plotting
    z_matrix = obj_df.pivot(index='y', columns='x', values='z').values
    x_unique = np.sort(obj_df['x'].unique())
    y_unique = np.sort(obj_df['y'].unique())

    return z_matrix, x_unique, y_unique

def plot_acqf_surface(obj_df, xaxis, yaxis, **kwargs):
    # Pivot the DataFrame to create matrices for plotting
    z_matrix, x_unique, y_unique = transform_obj_for_plot(obj_df)

    # Create the 3D surface plot with hover text
    surface = go.Figure(data=[go.Surface(
        z=z_matrix,
        x=x_unique,
        y=y_unique,
        hovertemplate=f"<b>{xaxis}</b>: %{{x}}<br><b>{yaxis}</b>: %{{y}}<br><b>Acquisition Function</b>: %{{z}}<extra></extra>"
    )])

    surface.update_layout(
        scene=dict(
            xaxis_title=xaxis,
            yaxis_title=yaxis,
            zaxis_title='Acquisition Function',
        ),
        width=800,
        height=600,
        **kwargs
    )
    surface.show()

def plot_acqf_contour(obj_df, xaxis, yaxis, **kwargs):
    # Pivot the DataFrame to create matrices for plotting
    z_matrix, x_unique, y_unique = transform_obj_for_plot(obj_df)

    # Create the 2D contour plot with hover text
    contour = go.Figure(data=[go.Contour(
        z=z_matrix,
        x=x_unique,
        y=y_unique,
        hovertemplate=f"<b>{xaxis}</b>: %{{x}}<br><b>{yaxis}</b>: %{{y}}<br><b>Acquisition Function</b>: %{{z}}<extra></extra>"
    )])

    contour.update_layout(
        xaxis_title=xaxis,
        yaxis_title=yaxis,
        width=800,
        height=600,
        **kwargs
    )
    contour.show()
