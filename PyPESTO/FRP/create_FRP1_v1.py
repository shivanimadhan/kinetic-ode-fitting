from ast import mod
import numpy as np
import pandas as pd
import petab
from petab.v1.C import PARAMETER_ID, PARAMETER_SCALE, LOWER_BOUND, UPPER_BOUND, NOMINAL_VALUE, ESTIMATE, OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT, OBSERVABLE_FORMULA, NOISE_FORMULA, LOG, CONDITION_ID, CONDITION_NAME, FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES
import amici
from libsbml import Model
import pypesto

import PyPESTO.FRP.sbml as sbml
import os

# from Sensitivity.FRP_Model_ import CopolymerizationModel
from PyPESTO.FRP.FRP_Model import CopolymerizationModel

MODEL_NAME = 'FRP1_v1'

kd_true = 1e-3
kp_true = 1e4
kdp_true = 1e2
kt_true = 1e6

kp_kt_ratio_true = kp_true / kt_true
kd_kt_true = kd_true * kt_true
Kr_true = kdp_true / kp_true


def get_rate_params_from_ratios(kpAA, rA, rB, rX, KBB, kt):
    
    kpAB = kpAA / rA
    kpBB = kpAA / rX
    kpBA = kpBB / rB
    
    kdBB = kpBB*KBB
    
    return kpAA, kpAB, kpBA, kpBB, kdBB, kt

def get_ratios_from_rate_params(kpAA, kpAB, kpBA, kpBB, kdBB, kt):
    
    rA = kpAA / kpAB
    rB = kpBB / kpBA
    rX = kpAA / kpBB
    KBB = kdBB / kpBB
    
    return kpAA, rA, rB, rX, KBB, kt

def define_parameters():
    # parameters_df = pd.DataFrame(
    #     data={
    #         PARAMETER_ID: ["rA", "rB", "rX", "kpAA", "KBB", "kt"],
    #         PARAMETER_SCALE: [LOG] * 6,
    #         LOWER_BOUND: [1e-3] * 3 + [1e3] + [1e-3] + [1e4] * 1,
    #         UPPER_BOUND: [1e3] * 3 + [1e6] + [1e3] + [1e8] * 1,
    #         NOMINAL_VALUE: [1.0, 1.0, 1.0, 1.0, 1.0, 1e6],
    #         ESTIMATE: [1, 1, 1, 1, 0, 1],
    #     }
    # ).set_index(PARAMETER_ID)
    # parameters_df = pd.DataFrame(
    #     data={
    #         PARAMETER_ID: ['kd', 'kp', 'kt'],
    #         PARAMETER_SCALE: [LOG] * 3,
    #         LOWER_BOUND: [1e-6, 1e3, 1e3],
    #         UPPER_BOUND: [1e3, 1e6, 1e8],
    #         NOMINAL_VALUE: [1e-4, 0, 0],
    #         ESTIMATE: [0, 1, 1],
    #     }
    # ).set_index(PARAMETER_ID)
    # parameters_df = pd.DataFrame(
    #     data={
    #         PARAMETER_ID: ['kp', 'kp_kt_ratio', 'kd_kt'],
    #         PARAMETER_SCALE: [LOG] * 3,
    #         LOWER_BOUND: [1e3, 1e-6, 1e-6],
    #         UPPER_BOUND: [1e8, 1e6, 1e6],
    #         NOMINAL_VALUE: [1e-4, 0, 0],
    #         ESTIMATE: [1, 1, 1],
    #     }
    # ).set_index(PARAMETER_ID)
    # return parameters_df
    parameters_df = pd.DataFrame(
        data={
            PARAMETER_ID: ['kp', 'kp_kt_ratio', 'kd_kt', 'Kr'],
            PARAMETER_SCALE: [LOG] * 4,
            LOWER_BOUND: [1e3, 1e-6, 1e-6, 1e-3],
            UPPER_BOUND: [1e8, 1e6, 1e6, 1e1],
            NOMINAL_VALUE: [1e-4, 0, 0, 0],
            ESTIMATE: [1, 1, 1, 1],
        }
    ).set_index(PARAMETER_ID)
    return parameters_df

def load_sbml_model(sbml_model_filepath):
    return sbml.create_SBML_FRP1_v1(
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
    constant_parameters = ["f"]
    observables = {
        "obs_m": {'formula': 'M'}}
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

def run_amici_simulation(model, timepoints, cI0: float=5e-3, cM0:float=0.0, sigma:float=0.0):
    
    # model
    print(model.getParameters())
    print(model.getParameterNames())
    
    # model.setParameterByName('kd', kd_true)
    model.setParameterByName('kp', kp_true)
    model.setParameterByName('kp_kt_ratio', kp_kt_ratio_true)
    model.setParameterByName('kd_kt', kd_kt_true)
    model.setParameterByName('Kr', Kr_true)
    # model.setParameterByName('kt', kt_true)
    
    # model.setParameterByName('kpAA', kpAA_true)
    # model.setParameterByName('rA', rA_true)
    # model.setParameterByName('rB', rB_true)
    # model.setParameterByName('rX', rX_true)
    # model.setParameterByName('KBB', KBB_true)
    
    # model.setParameterByName('kt', kt_true)
    # model.setParameterByName('ktAA', ktAA_true)
    # model.setParameterByName('ktBB', ktBB_true)
    # model.setParameterByName('ktAB', ktAB_true)

    print(model.getParameters())

    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-10)
    
    # timepoints = np.linspace(0, 100, 10)
    model.setInitialStates([cI0, 0, cM0, 0, 0])
    
    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)
    
    meas_m = rdata.by_id('M') + sigma * np.random.randn(len(timepoints))
    
    return meas_m
    
# def run_ode_simulation(timepoints, cI0: float=5e-3, cA0:float=0.0, cB0:float=0.0, sigma:float=0.0):    
    
#     # ktAB_true = np.sqrt(ktAA_true * ktBB_true)
    
#     k = [
#         kd_true, # Initiator dissociation rate constant
#         kpAA_true, kpAB_true, kpBA_true, kpBB_true, # Propagation rate constants
#         # params[0], params[1], params[2], params[3], # Propagation rate constants
#         0, 0, 0, kdBB_true, # Depropagation rate constants
#         0, kt_true,
#         0, kt_true,
#         0, kt_true,
#         # 0, 0, 0, # Termination (combination) rate constants
#         # ktAA_true, ktAB_true, ktBB_true  # Termination (disproportionation) rate constants
#     ]

#     y0 = np.zeros(33)
#     y0[0] = cI0
#     y0[2] = cA0
#     y0[3] = cB0

#     # Initialize and solve the model
#     cm = CopolymerizationModel(k=k, f=0.5, y0=y0)
    
#     # t0 = min(timepoints)
#     t0 = 0
#     t1 = max(timepoints)
    
#     # cm.solve([t0, t1], num_points=len(timepoints))
#     cm.solve([t0, t1], num_points=len(timepoints), remove_zero=False, solve_chain_model=False, solve_sequence_model=False)
    
#     # output_names = ['[A]', '[B]', 'nAvgCL', 'nAvgSL A', 'nAvgSL B']
#     # output_matrix[i] = np.array([cm.cA, cm.cB, cm.NACL, cm.NASL_A, cm.NASL_B]).T
    
#     # return t, output_matrix, output_names
#     return cm.cA, cm.cB
    
def define_FRP_measurements(amici_model):
    
    timepoints = np.linspace(0, (50*60), 20)
    
    obs_sigma = 0.02
    meas_sigma = 0.02
    
    # num_conditions = 6
    # # cI0  = [0.02, 0.02, 0.02]
    # cI0s = [0.005, 0.005, 0.01, 0.01, 0.02, 0.02]
    # # fA0s = [0.25, 0.5, 0.75, 0.25, 0.5, 0.75]
    # cM0s = [2.0, 1.0, 2.0, 1.0, 2.0, 1.0]
    
    num_conditions = 2
    cI0s = [0.01, 0.02]
    cM0s = [2.0, 2.0]
    
    # num_conditions = 3
    # # cI0  = [0.02, 0.02, 0.02]
    # cI0 = [0.005, 0.005, 0.005]
    # fA0s = [0.25, 0.5, 0.75]
    # cM0s = [2.0, 2.0, 2.0]
    
    # num_conditions = 2
    # fA0s = [0.75, 0.25]
    # cM0s = [3.0, 3.0]
     
    assert(len(cI0s) == len(cM0s) and len(cM0s) == num_conditions)
    
    condition_ids = [f'c_{i}' for i in range(num_conditions)]
    conditions = [(cI0s[i], cM0s[i]) 
        for i in range(num_conditions)
    ]
    
    obervables_df = pd.DataFrame(
        data={
            OBSERVABLE_ID: ['obs_m'],
            OBSERVABLE_FORMULA: ['M'],
            NOISE_FORMULA: [obs_sigma]
        }
    ).set_index(OBSERVABLE_ID)
    
    measurement_df = pd.DataFrame(columns=[OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT])
    conditions_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'I', 'M'])
    
    for i, (cI0, cM0) in enumerate(conditions):
        # print(i, cI0, cM0)
        
        meas_m = run_amici_simulation(amici_model, timepoints, cI0=cI0, cM0=cM0, sigma=meas_sigma)
        # print(meas_m)
        
        # meas_a_ode, meas_b_ode = run_ode_simulation(timepoints, cI0=cI0, cA0=cA0, cB0=cB0, sigma=meas_sigma)
        # print(meas_a_ode, meas_b_ode)
    
    
        df_m = pd.DataFrame({
            OBSERVABLE_ID: ['obs_m'] * len(timepoints),
            SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
            TIME: timepoints,
            MEASUREMENT: meas_m
        })
        
        conditions_data = pd.DataFrame({
            CONDITION_ID: [condition_ids[i]],
            CONDITION_NAME: [condition_ids[i]],
            'I': [cI0],
            'M': [cM0],
        })
        
        if i == 0:
            measurement_df = pd.concat([df_m], ignore_index=True)
            conditions_df = pd.concat([conditions_data], ignore_index=True)
        else:
            measurement_df = pd.concat([measurement_df, df_m], ignore_index=True)
            conditions_df = pd.concat([conditions_df, conditions_data], ignore_index=True)
        
        
    return obervables_df, conditions_df, measurement_df

# problem --------------------------------------------------------------------

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