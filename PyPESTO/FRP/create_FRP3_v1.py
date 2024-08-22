import numpy as np
import pandas as pd
import petab
from petab.v1.C import PARAMETER_ID, PARAMETER_SCALE, LOWER_BOUND, UPPER_BOUND, NOMINAL_VALUE, ESTIMATE, OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT, OBSERVABLE_FORMULA, NOISE_FORMULA, LOG, CONDITION_ID, CONDITION_NAME, FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES
import amici
from libsbml import Model

import PyPESTO.FRP.sbml as sbml
import os

MODEL_NAME = 'FRP3_v1'

kpAA_true = 2.0
kpAB_true = 0.5
kpAC_true = 0.25

kpBA_true = 10.0
kpBB_true = 5.0
kpBC_true = 2.5

kpCA_true = 5.0
kpCB_true = 2.5
kpCC_true = 1.25

rAB_true = kpAA_true / kpAB_true
rAC_true = kpAA_true / kpAC_true
rBA_true = kpBB_true / kpBA_true
rBC_true = kpBB_true / kpBC_true
rCA_true = kpCC_true / kpCA_true
rCB_true = kpCC_true / kpCB_true
rAxB_true = kpAA_true / kpBB_true
rAxC_true = kpAA_true / kpCC_true

# rAA_true = kpAA_true / kpAB_true
# rB_true = kpBB_true / kpBA_true
# rX_true = kpAB_true / kpAA_true


def define_parameters():
    
    parameters_df = pd.DataFrame(
        data={
            PARAMETER_ID: ["rAB", "rAC", "rBA", "rBC", "rCA", "rCB", "rAxB", "rAxC", "kpAA"],
            PARAMETER_SCALE: [LOG] * 9,
            LOWER_BOUND: [1e-3] * 9,
            UPPER_BOUND: [1e3] * 9,
            NOMINAL_VALUE: [1.0] * 9,
            ESTIMATE: [1] * 9,
        }
    ).set_index(PARAMETER_ID)
    return parameters_df
    


def load_sbml_model(sbml_model_filepath):
    return sbml.create_SBML_FRP3_v1(
        sbml_model_filepath,
        with_rules=False,
    )
    
def load_amici_from_sbml():
    
    sbml_model_filepath = f'/SBML/PyPESTO/FRP/{MODEL_NAME}/sbml_model.xml'
    model_dir = os.path.dirname(sbml_model_filepath)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    
    sbml_model = load_sbml_model(sbml_model_filepath)
    
    print('Importing AMICI model from SBML')
    sbml_importer = amici.SbmlImporter(sbml_model_filepath)
    
    model_output_dir = os.path.join('tmp/', MODEL_NAME)
    constant_parameters = ["kd", "f"]
    observables = {"obs_a": {'formula': 'A'},
                "obs_b": {'formula': 'B'},
                "obs_c": {'formula': 'C'}}
    sbml_importer.sbml2amici(
        MODEL_NAME, 
        model_output_dir, 
        verbose=False,
        observables=observables,
        constant_parameters=constant_parameters
    )
    
    model_module = amici.import_model_module(MODEL_NAME, model_output_dir)

    # Instantiate model
    model = model_module.getModel()
    
    return model, sbml_model_filepath

def run_amici_simulation(model, timepoints, cI0: float=5e-3, cA0:float=0.0, cB0:float=0.0, cC0:float=0.0, sigma:float=0.0):
    
    # model
    model.setParameterByName('kpAA', kpAA_true)
    model.setParameterByName('rAB', rAB_true)
    model.setParameterByName('rAC', rAC_true)
    model.setParameterByName('rBA', rBA_true)
    model.setParameterByName('rBC', rBC_true)
    model.setParameterByName('rCA', rCA_true)
    model.setParameterByName('rCB', rCB_true)
    model.setParameterByName('rAxB', rAxB_true)
    model.setParameterByName('rAxC', rAxC_true)
    # model.setParameterByName('rA', rA_true)
    # model.setParameterByName('rB', rB_true)
    # model.setParameterByName('rX', rX_true)

    # model.setParameterByName('kpBA', kpBA_true)
    # model.setParameterByName('kpBB', kpBB_true)
    print(model.getParameters())

    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-10)
    
    # timepoints = np.linspace(0, 100, 10)
    model.setInitialStates([cI0, 0, cA0, cB0, cC0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    
    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)
    
    meas_a = rdata.by_id('A') + sigma * np.random.randn(len(timepoints))
    meas_b = rdata.by_id('B') + sigma * np.random.randn(len(timepoints))
    meas_c = rdata.by_id('C') + sigma * np.random.randn(len(timepoints))
    
    return meas_a, meas_b, meas_c
    

    
def define_FRP_measurements(amici_model):
    
    timepoints = np.linspace(0, 100, 10)
    
    obs_sigma = 0.02
    meas_sigma = 0.00
    
    num_conditions = 6
    fA0s = [0.25, 0.5, 0.75, 0.40, 0.25, 0.15]
    fB0s = [0.45, 0.40, 0.10, 0.25, 0.5, 0.75]
    cM0s = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
    
    # num_conditions = 2
    # fA0s = [0.75, 0.25]
    # fB0s = [0.10, 0.20]
    # cM0s = [3.0, 3.0]
     
    assert(len(fA0s) == len(cM0s) and len(fA0s) == num_conditions)
    
    condition_ids = [f'c_{i}' for i in range(num_conditions)]
    conditions = [(fA0s[i]*cM0s[i], (fB0s[i])*cM0s[i], (1-fA0s[i]-fB0s[i])*cM0s[i]) 
        for i in range(num_conditions)
    ]
    
    obervables_df = pd.DataFrame(
        data={
            OBSERVABLE_ID: ['obs_a', 'obs_b', 'obs_c'],
            OBSERVABLE_FORMULA: ['A', 'B', 'C'],
            NOISE_FORMULA: [obs_sigma] * 3
        }
    ).set_index(OBSERVABLE_ID)
    
    measurement_df = pd.DataFrame(columns=[OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT])
    conditions_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'A', 'B', 'C'])
    
    for i, (cA0, cB0, cC0) in enumerate(conditions):
        meas_a, meas_b, meas_c = run_amici_simulation(amici_model, timepoints, cA0=cA0, cB0=cB0, cC0=cC0, sigma=meas_sigma)
    
        df_a = pd.DataFrame({
            OBSERVABLE_ID: ['obs_a'] * len(timepoints),
            SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
            TIME: timepoints,
            MEASUREMENT: meas_a
        })
        
        df_b = pd.DataFrame({
            OBSERVABLE_ID: ['obs_b'] * len(timepoints),
            SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
            TIME: timepoints,
            MEASUREMENT: meas_b
        })
        
        df_c = pd.DataFrame({
            OBSERVABLE_ID: ['obs_c'] * len(timepoints),
            SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
            TIME: timepoints,
            MEASUREMENT: meas_c
        })
        
        measurement_df = pd.concat([measurement_df, df_a, df_b, df_c], ignore_index=True)
        
        conditions_data = pd.DataFrame({
            CONDITION_ID: [condition_ids[i]],
            CONDITION_NAME: [condition_ids[i]],
            'A': [cA0],
            'B': [cB0],
            'C': [cC0],
        })
        
        conditions_df = pd.concat([conditions_df, conditions_data], ignore_index=True)
        
    return obervables_df, conditions_df, measurement_df

# problem --------------------------------------------------------------------

from petab.v2.lint import lint_problem
def write_petab_files(amici_model, sbml_model_filepath):
    
    parameters_df = define_parameters()
    observables_df, conditions_df, measurements_df = define_FRP_measurements(amici_model)
    
    model_dir = os.path.dirname(sbml_model_filepath)
    model_filename = os.path.basename(sbml_model_filepath)
        
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
    problem = petab.v2.Problem.from_yaml(yaml_filepath)
    print(lint_problem(problem))
    
    return yaml_filepath
