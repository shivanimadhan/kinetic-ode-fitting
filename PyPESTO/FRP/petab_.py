import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass

import os

from petab.v1.C import FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES

import logging

import numpy as np
import petab
import amici
import amici.plotting
import sys
import os
import importlib

import pypesto
import pypesto.optimize as optimize
import pypesto.petab
import pypesto.sample as sample
import pypesto.visualize as visualize

@dataclass
class PetabData:
    conditions: pd.DataFrame
    observables: pd.DataFrame
    measurements: pd.DataFrame
    parameters: pd.DataFrame



def plot_measurements(mdf: pd.DataFrame, conversion=False):
    # Get the unique conditions and observables
    conditions = mdf["simulationConditionId"].unique()
    observables = mdf["observableId"].unique()

    # Calculate the number of rows and columns for subplots
    num_conditions = len(conditions)
    num_cols = 3
    num_rows = (num_conditions + num_cols - 1) // num_cols  # Round up to get the number of rows needed

    # Create the subplots
    fig, ax = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    ax = ax.flatten()  # Flatten to easily iterate

    # Loop through each condition and plot the data
    for i, condition in enumerate(conditions):
        condition_data = mdf[mdf["simulationConditionId"] == condition]

        for obs in observables:
            obs_data = condition_data[condition_data['observableId'] == obs]
            
            if conversion:
                max_concentration = obs_data['measurement'].max()
                print('max', max_concentration)
                conv_data = (max_concentration - obs_data['measurement']) / max_concentration
                ax[i].plot(obs_data['time'], conv_data, 'o-', label=obs)
            else:
                ax[i].plot(obs_data['time'], obs_data['measurement'], 'o-', label=obs)

        ax[i].set_title(f'Condition {condition}')
        ax[i].set_xlabel('Time (s)')
        # ax[i].set_ylabel('Simulated Concentration (M)')
        ax[i].set_ylim([0, 2.5])
        ax[i].legend(frameon=False)
        
        if conversion:
            ax[i].set_ylabel('Conversion')
            ax[i].set_ylim([-0.1, 1.1])
        else:
            ax[i].set_ylabel('Simulated Concentration (M)')
            ax[i].set_ylim([0, mdf['measurement'].max()])

    # Hide any unused subplots
    for j in range(i + 1, len(ax)):
        fig.delaxes(ax[j])

    plt.tight_layout()
    return fig, ax

def load_petab_problem(yaml_filepath: str, model_name: str) -> pypesto.problem.Problem:
    
    logger = logging.getLogger("pypesto.sample.diagnostics")
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())

    # import to petab
    petab_problem = petab.Problem.from_yaml(
        yaml_filepath
    )

    importer = pypesto.petab.PetabImporter.from_yaml(
        yaml_filepath,
        model_name=model_name,
    )

    # Works when it runs twice for some reason
    try:
        problem = importer.create_problem(force_compile=True)
    except:
        problem = importer.create_problem(force_compile=True)
        
    return problem



def write_petab_files(petab_data: PetabData, sbml_model_filepath: str, model_name: str = ''):
    
    model_dir = os.path.dirname(sbml_model_filepath)
    model_filename = os.path.basename(sbml_model_filepath)
        
    petab.v1.write_condition_df(petab_data.conditions,     os.path.join(model_dir, "conditions.tsv"))
    petab.v1.write_measurement_df(petab_data.measurements, os.path.join(model_dir, "measurements.tsv"))
    petab.v1.write_observable_df(petab_data.observables,   os.path.join(model_dir, "observables.tsv"))
    petab.v1.write_parameter_df(petab_data.parameters,     os.path.join(model_dir, "parameters.tsv"))
    
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
    
    yaml_filepath = os.path.join(model_dir, f'{model_name}.yaml')
    petab.v1.yaml.write_yaml(yaml_config, yaml_filepath)

    # validate written PEtab files
    problem = petab.v1.Problem.from_yaml(yaml_filepath)
    petab.v1.lint.lint_problem(problem)
    
    return yaml_filepath