from SALib.analyze.sobol import analyze as sobol_analyze
from SALib.sample.sobol import sample as sobol_sample

import numpy as np
from tqdm import tqdm
from FRP_Model import CopolymerizationModel

from matplotlib import pyplot as plt
import pandas as pd
import os
import json

import warnings
# Suppress all warnings
warnings.filterwarnings("ignore")

### Sampling, solving, and analyzing for sensitivity analysis ###

def log_sample(problem, num_samples=1000) -> np.ndarray:
    
    n = len(problem['bounds'])
    log_bounds = []
    
    for i in range(n-1):
        log_bounds.append(np.log(problem['bounds'][i]))

    log_bounds.append(problem['bounds'][n-1])

    #log_bounds = [list(np.log(bounds)) for bounds in problem['bounds']]
    
    log_problem = {
        'num_vars': problem['num_vars'],
        'names': problem['names'],
        'bounds': log_bounds
    }
    
    log_param_values = sobol_sample(log_problem, num_samples)
    
    param_values = log_param_values.copy()
    param_values[:, 0:n-1] = np.exp(param_values[:, 0:n-1])

    return param_values

def solve_odes(param_values, time_points=20, tmax=10) -> np.ndarray:
    
    t = np.linspace(0, tmax, time_points)
    
    num_input_params = 5
    num_samples = param_values.shape[0]

    output_matrix = np.empty([num_samples, time_points, num_input_params])
    
    for i, params in tqdm(enumerate(param_values)):
        
        k = [
            1e8, 0.5, # Initiator dissociation rate constant
            params[0], params[1], params[2], params[3], # Propagation rate constants
            0, 0, 0, params[4], # Depropagation rate constants
            0, 0, 0, # Termination (combination) rate constants
            0, 0, 0  # Termination (disproportionation) rate constants
        ]
        
        fA = params[5]
        M0 = 3.0
        
        # Define initial concentrations
        I0 = 0.005
        A0 = fA * M0
        B0 = (1 - fA) * M0

        y0 = np.zeros(33)
        y0[0] = I0
        y0[2] = A0
        y0[3] = B0
    
        # Initialize and solve the model
        cm = CopolymerizationModel(k, y0)
        cm.solve([0, tmax], num_points=time_points)
        
        output_names = ['[A]', '[B]', 'nAvgCL', 'nAvgSL A', 'nAvgSL B']
        output_matrix[i] = np.array([cm.cA, cm.cB, cm.NACL, cm.NASL_A, cm.NASL_B]).T
    
    return t, output_matrix, output_names

def add_data_to_large_df(time_point, small_df, large_df):
    # Initialize a dictionary to store data to be appended
    data = {'time_point': time_point}
    
    for idx in small_df.index:
        data[f'{idx}_ST'] = small_df.loc[idx, 'ST']
        data[f'{idx}_ST_conf'] = small_df.loc[idx, 'ST_conf']

    df_row = pd.DataFrame([data])
    
    # Filter out empty or all-NA rows before concatenation
    if not df_row.empty and not df_row.isna().all(axis=1).all():
        large_df = pd.concat([large_df, df_row], ignore_index=True)
    
    return large_df

def sensitivity_analysis(problem, t, output_matrix, output_dir):
    
    num_outputs = output_matrix.shape[2]
    dfs = []
    
    for output in range(num_outputs):
        print('Performing sensitivity analysis on output', output, '...')
        
         # Create large dataframe
        df_large = pd.DataFrame(columns=['time_point'])
        
        for i, t_point in tqdm(enumerate(t)):
            if (i == 0):
                continue

            Si = sobol_analyze(problem, output_matrix[:, i, output])
            
            ST_df = Si.to_df()[0]

            df_large = add_data_to_large_df(t_point, ST_df, df_large)
            
        output_filepath = os.path.join(output_dir, f'{output}.csv')
        df_large.to_csv(output_filepath, index = False)
        dfs.append(df_large)
            
    return dfs

### Loading and saving data ###

def create_output_dir(output_dir: str):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Directory '{output_dir}' created.")
    else:
        print(f"Directory '{output_dir}' already exists.")

def save_problem(problem, output_dir: str):
    output_filepath = os.path.join(output_dir, 'problem.json')
    with open(output_filepath, 'w') as f:
        json.dump(problem, f)

def save_params(param_values, output_dir: str):
    output_filepath = os.path.join(output_dir, 'params.npy')
    np.save(output_filepath, param_values)
    
def save_ode_outputs(output_mat, t, output_dir: str):
    output_filepath = os.path.join(output_dir, 'outputs.npy')
    np.save(output_filepath, output_mat)
    
    output_filepath = os.path.join(output_dir, 'times.npy')
    np.save(output_filepath, t)
    
def is_csv_file(filepath):
    if not os.path.exists(filepath):
        return False
    
    if not os.path.isfile(filepath):
        return False
    
    if os.path.splitext(filepath)[1].lower() != '.csv':
        return False
    
    return True

def load_data(input_dir):

    csv_filepaths = [os.path.join(input_dir, name) for name in os.listdir(input_dir) if is_csv_file(os.path.join(input_dir, name))]
    csv_filepaths = sorted(csv_filepaths)

    # Printing the file paths
    dfs = []
    for filepath in csv_filepaths:
        df = pd.read_csv(filepath)
        dfs.append(df)
        
    return dfs
   
def load_params(dir: str):
    filepath = os.path.join(dir, 'params.npy')
    return np.load(filepath)

def load_ode_outputs(dir: str):
    filepath = os.path.join(dir, 'outputs.npy')
    output_mat = np.load(filepath)
    
    filepath = os.path.join(dir, 'times.npy')
    t = np.load(filepath)
    return t, output_mat

def load_problem(dir: str):
    filepath = os.path.join(dir, 'problem.json')
    with open(filepath, 'r') as f:
        problem = json.load(f)
    return problem
    
### Plotting functions ###

def plot_ode_outputs(output_mat: np.ndarray, output_names: list, t=None):
    
    num_sets = output_mat.shape[0]

    fig, axs = plt.subplots(1, 5, figsize=(25, 5))
    
    for i in range(len(axs)):
        for p in range(num_sets):
            if t is None:
                axs[i].plot(output_mat[p, :, i], '-', alpha=0.05)
            else:   
                axs[i].plot(t, output_mat[p, :, i], '-', alpha=0.05)
            
        axs[i].tick_params(labelsize=18)
        axs[i].set_xlabel('Time (s)', fontsize=18)
        axs[i].set_ylabel(output_names[i], fontsize=18)

    # Format axes
    axs[0].set_yticks([0, 0.5, 1.0, 1.5])
    axs[0].set_yticklabels([0, 0.5, 1.0, 1.5])
    
    axs[1].set_yticks([0, 0.5, 1.0, 1.5])
    axs[1].set_yticklabels([0, 0.5, 1.0, 1.5])
    
    axs[2].set_yticks([0, 200, 400, 600])
    axs[2].set_yticklabels([0, 200, 400, 600])


    axs[3].set_yscale('log')
    axs[3].set_yticks([1, 2, 4, 8, 16, 32, 64])
    axs[3].set_yticklabels([1, 2, 4, 8, 16, 32, 64])
    
    axs[4].set_yscale('log')
    axs[4].set_yticks([1, 2, 4, 8, 16, 32, 64])
    axs[4].set_yticklabels([1, 2, 4, 8, 16, 32, 64])
    
    plt.tight_layout()
    
    return fig, axs

def plot_param_values(param_values, idx1=0, idx2=1, names=None):
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    ax.scatter(param_values[:,idx1], param_values[:,idx2], s=0.1, alpha=0.5)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    if names:
        ax.set_xlabel(names[idx1], fontsize=14)
        ax.set_ylabel(names[idx2], fontsize=14, rotation=90)

    ax.tick_params(labelsize=12)

    ax.set_title('Sampled Rate Constants', fontsize=16)
    
    plt.tight_layout()
    
    return fig, ax

def plot_sensitivity(ax, df, param_number, data_name, conf_name=None, color='k'):
    
    # Plot the line graph
    ax.plot(df['time_point'], df[data_name], label=data_name, color=color, linewidth=4)

    # Shade the confidence interval
    if conf_name is not None:
        ax.fill_between(df['time_point'], df[data_name] - df[conf_name], 
            df[data_name] + df[conf_name], color=color, alpha=0.3)
    
    # Show plot
    ax.grid(True)
    
    return ax

def plot_sensitivity_grid(problem, dfs):
    
    num_outputs = len(dfs)
    num_vars = problem['num_vars']
    
    print(num_outputs)
    print(num_vars)

    sz = 6
    fig, axs = plt.subplots(num_outputs, num_vars, figsize=(num_vars*sz, num_outputs*sz), dpi=300)

    dtype = 'ST'
    data_labels = [name+'_'+dtype for name in problem['names']]
    conf_labels = [name+'_'+dtype+'_conf' for name in problem['names']]
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
    output_names = ['[A]', '[B]', 'nAvgCL', 'nAvgSL A', 'nAvgSL B']

    for i, df in enumerate(dfs): # Loop through each output variable
        for j, data_label in enumerate(data_labels): # Loop through each parameter
            c = colors[i]
            n = output_names[i]

            ax = axs[i, j]
            
            ax = plot_sensitivity(ax, df, j, data_label, conf_name=conf_labels[j], color=c)
            ax.set_ylim([0, 1])
            ax.legend([n], loc='upper right', fontsize=36)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
            
            ax.set_xticks([0, 50, 100, 150, 200])
            ax.tick_params(labelsize=36)
            
            if i == 0:
                ax.set_title(f'{problem["names"][j]}', fontsize=48, fontweight='bold')
                
            if j == 0:
                ax.set_ylabel(output_names[i], fontsize=48, fontweight='bold')
                ax.set_yticklabels([0, '', 0.5, '', 1])
            else:
                ax.set_yticklabels([])
                
            if i == num_outputs-1:
                ax.set_xlabel('Time (s)', fontsize=36)
                ax.set_xticklabels([0, '', 100, '', 200])
            else:
                ax.set_xticklabels([])
                
    plt.tight_layout()
    return fig, axs