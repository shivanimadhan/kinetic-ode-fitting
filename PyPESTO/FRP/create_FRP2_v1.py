import numpy as np
import pandas as pd
import petab
from petab.v1.C import PARAMETER_ID, PARAMETER_SCALE, LOWER_BOUND, UPPER_BOUND, NOMINAL_VALUE, ESTIMATE, OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT, OBSERVABLE_FORMULA, NOISE_FORMULA, LOG, CONDITION_ID, CONDITION_NAME, FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES
import amici
from libsbml import Model
import pypesto

import PyPESTO.FRP.sbml as sbml
import os

MODEL_NAME = 'FRP2_v1'

from PyPESTO.FRP.FRP_Model import CopolymerizationModel

kpAA_true = 2.0
kpAB_true = 0.5
kpBA_true = 10.0
kpBB_true = 5.0

rA_true = kpAA_true / kpAB_true
rB_true = kpBB_true / kpBA_true
rX_true = kpAA_true / kpBB_true
# rX_true = kpAB_true / kpAA_true

def get_rate_params_from_ratios(kpAA, rA, rB, rX):
    
    kpAB = kpAA / rA
    kpBB = kpAA / rX
    kpBA = kpBB / rB
    
    return kpAA, kpAB, kpBA, kpBB

def get_ratios_from_rate_params(kpAA, kpAB, kpBA, kpBB):
    
    rA = kpAA / kpAB
    rB = kpBB / kpBA
    rX = kpAA / kpBB
    
    return kpAA, rA, rB, rX


def define_parameters():
    
    parameters_df = pd.DataFrame(
        data={
            PARAMETER_ID: ["rA", "rB", "rX", "kpAA"],
            PARAMETER_SCALE: [LOG] * 4,
            LOWER_BOUND: [1e-3] * 4,
            UPPER_BOUND: [1e3] * 4,
            NOMINAL_VALUE: [1.0, 1.0, 1.0, 1.0],
            ESTIMATE: [1, 1, 1, 1],
        }
    ).set_index(PARAMETER_ID)
    return parameters_df

def load_sbml_model(sbml_model_filepath):
    return sbml.create_SBML_FRP2_v1(
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
                "obs_b": {'formula': 'B'},}
    sbml_importer.sbml2amici(
        MODEL_NAME, 
        model_output_dir, 
        verbose=False,
        observables=observables,
        constant_parameters=constant_parameters
    )
    print(f'Import model ({MODEL_NAME}) module.')
    model_module = amici.import_model_module(MODEL_NAME, model_output_dir)

    # Instantiate model
    model = model_module.getModel()
    
    return model, sbml_model_filepath

def run_amici_simulation(model, timepoints, cI0: float=5e-3, cA0:float=0.0, cB0:float=0.0, sigma:float=0.0):
    
    model.setParameterByName('kpAA', kpAA_true)
    model.setParameterByName('rA', rA_true)
    model.setParameterByName('rB', rB_true)
    model.setParameterByName('rX', rX_true)

    # print(model.getParameters())

    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-10)
    
    model.setInitialStates([cI0, 0, cA0, cB0, 0, 0, 0, 0, 0, 0])
    
    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)
    
    meas_a = rdata.by_id('A') + sigma * np.random.randn(len(timepoints))
    meas_b = rdata.by_id('B') + sigma * np.random.randn(len(timepoints))
    # print(rdata.by_id('I'))
    
    return meas_a, meas_b
    
def run_ode_simulation(timepoints, cI0: float=5e-3, cA0:float=0.0, cB0:float=0.0, sigma:float=0.0):    

    k = [
        1e10, # Initiator dissociation rate constant
        kpAA_true, kpAB_true, kpBA_true, kpBB_true, # Propagation rate constants
        # params[0], params[1], params[2], params[3], # Propagation rate constants
        0, 0, 0, 0, # Depropagation rate constants
        0, 0, 
        0, 0,
        0, 0,
    ]
    
    # print(kpAA_true, kpAB_true, kpBA_true, kpBB_true)

    y0 = np.zeros(33)
    y0[0] = cI0
    y0[2] = cA0
    y0[3] = cB0

    # Initialize and solve the model
    cm = CopolymerizationModel(k=k, f=0.5, y0=y0)
    
    t0 = 0
    t1 = max(timepoints)
    
    cm.solve([t0, t1], num_points=len(timepoints), remove_zero=False, solve_chain_model=False, solve_sequence_model=False)
    
    # output_names = ['[A]', '[B]', 'nAvgCL', 'nAvgSL A', 'nAvgSL B']
    # output_matrix[i] = np.array([cm.cA, cm.cB, cm.NACL, cm.NASL_A, cm.NASL_B]).T
    
    # return t, output_matrix, output_names
    # print(cm.cI)
    return cm.t, cm.cA, cm.cB
       

def define_FRP_measurements(amici_model):
    
    timepoints = np.linspace(0, 100, 10)
    
    # num_conditions = 2
    # # fA0s = [0.25, 0.5]
    # fA0s = [0.5, 0.75]
    # cM0s = [3.0, 3.0]
    
    num_conditions = 3
    fA0s = [0.25, 0.5, 0.75]
    cI0s = [0.005, 0.005, 0.005]
    cM0s = [3.0, 3.0, 3.0]
    
    obs_sigma = 0.02
    meas_sigma = 0.00
     
    assert(len(fA0s) == len(cM0s) and len(fA0s) == num_conditions)
    
    condition_ids = [f'c_{i}' for i in range(num_conditions)]
    conditions = [(cI0s[i], fA0s[i]*cM0s[i], (1-fA0s[i])*cM0s[i]) 
        for i in range(num_conditions)
    ]
    
    obervables_df = pd.DataFrame(
        data={
            OBSERVABLE_ID: ['obs_a', 'obs_b'],
            OBSERVABLE_FORMULA: ['A', 'B'],
            NOISE_FORMULA: [obs_sigma] * 2
        }
    ).set_index(OBSERVABLE_ID)
    
    measurement_df = pd.DataFrame(columns=[OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT])
    conditions_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'R', 'A', 'B'])
    
    for i, (cI0, cA0, cB0) in enumerate(conditions):
        meas_a, meas_b = run_amici_simulation(amici_model, timepoints, cI0=cI0, cA0=cA0, cB0=cB0, sigma=meas_sigma)
        t_ode, meas_a_ode, meas_b_ode = run_ode_simulation(timepoints, cI0=cI0, cA0=cA0, cB0=cB0, sigma=meas_sigma)
        
        tol = 1e-8
        # assert(np.allclose(meas_a, meas_a_ode, atol=tol))
        # assert(np.allclose(meas_b, meas_b_ode, atol=tol))
        
        # print(f'SBML {i}:', meas_a, meas_b)
        # print(f'ODE  {i}:', meas_a_ode, meas_b_ode)
    
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

        conditions_data = pd.DataFrame({
            CONDITION_ID: [condition_ids[i]],
            CONDITION_NAME: [condition_ids[i]],
            'R': [cI0],
            'A': [cA0],
            'B': [cB0]
        })
        
        if i == 0:
            measurement_df = pd.concat([df_a, df_b], ignore_index=True)
            conditions_df = pd.concat([conditions_data], ignore_index=True)
        else:
            measurement_df = pd.concat([measurement_df, df_a, df_b], ignore_index=True)
            conditions_df = pd.concat([conditions_df, conditions_data], ignore_index=True)
        
    return obervables_df, conditions_df, measurement_df

# problem --------------------------------------------------------------------

import matplotlib.pyplot as plt



# def plot_measurements(mdf: pd.DataFrame):
    

#     fig, ax = plt.subplots(2, 3, figsize=(15, 10))


#     c0_data = mdf[mdf["simulationConditionId"] == "c_0"]
#     c1_data = mdf[mdf["simulationConditionId"] == "c_1"]
#     c2_data = mdf[mdf["simulationConditionId"] == "c_2"]
#     c3_data = mdf[mdf["simulationConditionId"] == "c_3"]
#     c4_data = mdf[mdf["simulationConditionId"] == "c_4"]
#     c5_data = mdf[mdf["simulationConditionId"] == "c_5"]
#     # print(c0_data)

#     ax[0][0].plot(c0_data[c0_data['observableId']=='obs_a']['time'], c0_data[c0_data['observableId']=='obs_a']['measurement'], 'o-', label='obs_a')
#     ax[0][0].plot(c0_data[c0_data['observableId']=='obs_b']['time'], c0_data[c0_data['observableId']=='obs_b']['measurement'], 'o-', label='obs_b')
#     ax[0][0].plot(c0_data[c0_data['observableId']=='obs_c']['time'], c0_data[c0_data['observableId']=='obs_c']['measurement'], 'o-', label='obs_c')

#     ax[0][1].plot(c1_data[c1_data['observableId']=='obs_a']['time'], c1_data[c1_data['observableId']=='obs_a']['measurement'], 'o-', label='obs_a')
#     ax[0][1].plot(c1_data[c1_data['observableId']=='obs_b']['time'], c1_data[c1_data['observableId']=='obs_b']['measurement'], 'o-', label='obs_b')
#     ax[0][1].plot(c1_data[c1_data['observableId']=='obs_c']['time'], c1_data[c1_data['observableId']=='obs_c']['measurement'], 'o-', label='obs_c')

#     ax[0][2].plot(c2_data[c2_data['observableId']=='obs_a']['time'], c2_data[c2_data['observableId']=='obs_a']['measurement'], 'o-', label='obs_a')
#     ax[0][2].plot(c2_data[c2_data['observableId']=='obs_b']['time'], c2_data[c2_data['observableId']=='obs_b']['measurement'], 'o-', label='obs_b')
#     ax[0][2].plot(c2_data[c2_data['observableId']=='obs_c']['time'], c2_data[c2_data['observableId']=='obs_c']['measurement'], 'o-', label='obs_c')

#     ax[1][0].plot(c3_data[c3_data['observableId']=='obs_a']['time'], c3_data[c3_data['observableId']=='obs_a']['measurement'], 'o-', label='obs_a')
#     ax[1][0].plot(c3_data[c3_data['observableId']=='obs_b']['time'], c3_data[c3_data['observableId']=='obs_b']['measurement'], 'o-', label='obs_b')
#     ax[1][0].plot(c3_data[c3_data['observableId']=='obs_c']['time'], c3_data[c3_data['observableId']=='obs_c']['measurement'], 'o-', label='obs_c')

#     ax[1][1].plot(c4_data[c4_data['observableId']=='obs_a']['time'], c4_data[c4_data['observableId']=='obs_a']['measurement'], 'o-', label='obs_a')
#     ax[1][1].plot(c4_data[c4_data['observableId']=='obs_b']['time'], c4_data[c4_data['observableId']=='obs_b']['measurement'], 'o-', label='obs_b')
#     ax[1][1].plot(c4_data[c4_data['observableId']=='obs_c']['time'], c4_data[c4_data['observableId']=='obs_c']['measurement'], 'o-', label='obs_c')

#     ax[1][2].plot(c5_data[c5_data['observableId']=='obs_a']['time'], c5_data[c5_data['observableId']=='obs_a']['measurement'], 'o-', label='obs_a')
#     ax[1][2].plot(c5_data[c5_data['observableId']=='obs_b']['time'], c5_data[c5_data['observableId']=='obs_b']['measurement'], 'o-', label='obs_b')
#     ax[1][2].plot(c5_data[c5_data['observableId']=='obs_c']['time'], c5_data[c5_data['observableId']=='obs_c']['measurement'], 'o-', label='obs_c')

#     ax[0][2].legend(frameon=False)

#     for a in ax.flat:
#         a.set(xlabel='time (s)', ylim=[0, 2.5])
#         a.legend(frameon=False)
#         # a.set()
        
#     ax[0][0].set_ylabel('Simulated Concentration (M)')
#     ax[1][0].set_ylabel('Simulated Concentration (M)')
    
#     return fig, ax

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