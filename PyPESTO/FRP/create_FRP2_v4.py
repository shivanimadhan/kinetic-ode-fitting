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

MODEL_NAME = 'FRP2_v4'

from dataclasses import dataclass

@dataclass
class FRP2_v4_Parameters:
    rA: float
    rB: float
    rX: float
    KAA: float
    KAB: float
    KBA: float
    KBB: float
    kpAA: float
    kt_kp_ratio: float
    kd_kt: float

rA_true = 0.5
rB_true = 18
rX_true = 1
KAA_true = 0.0
KAB_true = 0.0
KBA_true = 0.0
KBB_true = 0.5
kpAA_true = 5e4
kd_true = 1e-3
kt_true = 1e6

kt_kp_ratio_true = kt_true / kpAA_true
kd_kt_true = kd_true * kt_true
# kd_kp_true = kd_true * kpAA_true
kd_kp_ratio_true = kd_true / kpAA_true

def get_rate_params_from_ratios(kpAA, rA, rB, rX, KAA, KAB, KBA, KBB, kt):
    
    kpAB = kpAA / rA
    kpBB = kpAA / rX
    kpBA = kpBB / rB
    
    kdAA = kpAA*KAA
    kdAB = kpAB*KAB
    kdBA = kpBA*KBA
    kdBB = kpBB*KBB
    
    return kpAA, kpAB, kpBA, kpBB, kdAA, kdAB, kdBA, kdBB, kt

def get_ratios_from_rate_params(kpAA, kpAB, kpBA, kpBB, kdAA, kdAB, kdBA, kdBB, kt):
    
    rA = kpAA / kpAB
    rB = kpBB / kpBA
    rX = kpAA / kpBB
    KAA = kdAA / kpAA
    KAB = kdAB / kpAB
    KBA = kdBA / kpBA
    KBB = kdBB / kpBB
    
    return kpAA, rA, rB, rX, KAA, KAB, KBA, KBB, kt

# def define_parameters():
#     parameters_df = pd.DataFrame(
#         data={
#             PARAMETER_ID: ['rA', 'rB', 'rX', 'KAA', 'KAB', 'KBA', 'KBB', 'kpAA', 'kt_kp_ratio', 'kd_kp_ratio'],
#             PARAMETER_SCALE: [LOG] * 3 + [LIN] * 4 + [LOG] * 3,
#             LOWER_BOUND: [1e-2] * 3 + [0] * 4 + [1e0] + [1e-6] * 2,
#             UPPER_BOUND: [1e2] * 3  + [1] * 4 + [1e6] + [1e6] * 2,
#             NOMINAL_VALUE: [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
#             ESTIMATE: [1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
#             ESTIMATE: [1, 1, 0, 1, 1, 1, 1, 0, 0, 0],
#         }
#     ).set_index(PARAMETER_ID)
#     return parameters_df

def define_parameters():
    parameters_df = pd.DataFrame(
        data={
            PARAMETER_ID: ['rA', 'rB', 'rX', 'KAA', 'KAB', 'KBA', 'KBB', 'kpAA'],#, 'kt_kp_ratio', 'kd_kp_ratio'],
            PARAMETER_SCALE: [LOG10] * 3 + [LIN] * 4 + [LOG10] * 1,#[LOG] * 3,
            LOWER_BOUND: [1e-2] * 3 + [0] * 4 + [1e0], #+ [1e-6] * 2,
            UPPER_BOUND: [1e2] * 3  + [1] * 4 + [1e6], #+ [1e6] * 2,
            NOMINAL_VALUE: [1.0, 0.8, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],#, 0.0, 1.0],
            # ESTIMATE: [1, 1, 0, 0, 0, 0, 1, 0],#, 0, 0],
            # ESTIMATE: [1, 1, 0, 1, 1, 1, 1, 0],#, 0, 0],
            ESTIMATE: [1, 1, 0, 1, 0, 0, 0, 0]
        }
    ).set_index(PARAMETER_ID)
    return parameters_df
    

def load_sbml_model(sbml_model_filepath):
    return sbml.create_SBML_FRP2_v4(
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
        "obs_a": {'formula': 'A'},
        "obs_b": {'formula': 'B'},}
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

def run_amici_simulation(model, timepoints, cI0: float=5e-3, cA0:float=0.0, cB0:float=0.0, sigma:float=0.0):
    
    # model
    # print(model.getParameters())
    # print(model.getParameterNames())
    # model.setParameterByName('kpAA', kpAA_true)
    model.setParameterByName('kpAA', kpAA_true)
    model.setParameterByName('rA', rA_true)
    model.setParameterByName('rB', rB_true)
    model.setParameterByName('rX', rX_true)
    # model.setParameterByName('rX', 1)
    model.setParameterByName('KAA', KAA_true)
    model.setParameterByName('KAB', KAB_true)
    model.setParameterByName('KBA', KBA_true)
    model.setParameterByName('KBB', KBB_true)
    
    model.setParameterByName('kd', kd_true)
    model.setParameterByName('kt', kt_true)
    
    # model.setParameterByName('kt_kp_ratio', kt_kp_ratio_true)
    # model.setParameterByName('kd_kp_ratio', kd_kp_ratio_true)
    # model.setParameterByName('kd_kt', kd_kt_true)
    # model.setParameterByName('kd_kp', kd_kp_true)
    # model.setParameterByName('kt_kp_ratio', kt_true / kpAA_true)
    # model.setParameterByName('kd_kp', kd_true * kpAA_true)
    # print(model.getParameters())

    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-10)
    
    # timepoints = np.linspace(0, 100, 10)
    model.setInitialStates([0, cI0, cA0, cB0, 0, 0, 0, 0, 0, 0, 0])
    # model.setInitialStates([cI0, 0, cA0, cB0, 0, 0, 0, 0, 0, 0, 0])
    
    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)
    
    meas_a = rdata.by_id('A') + sigma * np.random.randn(len(timepoints))
    meas_b = rdata.by_id('B') + sigma * np.random.randn(len(timepoints))
    
    return meas_a, meas_b
    
def define_FRP_measurements(amici_model):
    
    # timepoints = np.linspace(0, (50*60), 20)
    timepoints = np.linspace(0, 500, 20)
    
    obs_sigma = 0.02
    meas_sigma = 0.005
    
    # num_conditions = 6
    # # cI0  = [0.02, 0.02, 0.02]
    # cI0 = [0.005, 0.005, 0.005, 0.02, 0.02, 0.02]
    # fA0s = [0.25, 0.5, 0.75, 0.25, 0.5, 0.75]
    # cM0s = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    
    # num_conditions = 3
    # # # cI0  = [0.02, 0.02, 0.02]
    # cI0 = [0.02, 0.02, 0.02]
    # fA0s = [0.10, 0.20, 0.30]
    # # fA0s = [0.90, 0.80, 0.70]
    # cM0s = [1.0, 1.0, 1.0]
    
    # num_conditions = 5
    # # # cI0  = [0.02, 0.02, 0.02]
    # cI0 = [0.02, 0.02, 0.02, 0.02, 0.02]
    # fA0s = [0.10, 0.25, 0.5, 0.75, 0.90]
    # # fA0s = [0.90, 0.80, 0.70]
    # cM0s = [1.0, 1.0, 1.0, 1.0, 1.0]
    
    num_conditions = 3
    # # cI0  = [0.02, 0.02, 0.02]
    cI0 = [0.02, 0.02, 0.02]
    fA0s = [0.25, 0.50, 0.75]
    # fA0s = [0.90, 0.80, 0.70]
    cM0s = [1.0, 1.0, 1.0]
    
    # num_conditions = 1
    # # # cI0  = [0.02, 0.02, 0.02]
    # cI0 = [0.02]
    # fA0s = [0.75]
    # # fA0s = [0.90, 0.80, 0.70]
    # cM0s = [1.0]
    
    # num_conditions = 2
    # fA0s = [0.75, 0.25]
    # cM0s = [3.0, 3.0]
     
    assert(len(fA0s) == len(cM0s) and len(fA0s) == num_conditions)
    
    condition_ids = [f'c_{i}' for i in range(num_conditions)]
    conditions = [(cI0[i], fA0s[i]*cM0s[i], (1-fA0s[i])*cM0s[i]) 
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
    conditions_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'I', 'A', 'B'])
    
    for i, (cI0, cA0, cB0) in enumerate(conditions):
        print(i, cI0, cA0, cB0)
        
        meas_a, meas_b = run_amici_simulation(amici_model, timepoints, cI0=cI0, cA0=cA0, cB0=cB0, sigma=meas_sigma)
        print(meas_a, meas_b)
        
        # meas_a_ode, meas_b_ode = run_ode_simulation(timepoints, cI0=cI0, cA0=cA0, cB0=cB0, sigma=meas_sigma)
        # print(meas_a_ode, meas_b_ode)
    
    
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


def generate_FRP_data(amici_model, num_points=100, real=False):
    
    timepoints = np.linspace(0, (50*60), num_points)
    
    obs_sigma = 0.02
    meas_sigma = 0.00
    
    num_conditions = 6
    # cI0  = [0.02, 0.02, 0.02]
    cI0 = [0.005, 0.005, 0.005, 0.02, 0.02, 0.02]
    fA0s = [0.25, 0.5, 0.75, 0.25, 0.5, 0.75]
    cM0s = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    
    # cI0 = [0.005, 0.005, 0.005]
    # fA0s = [0.10, 0.20, 0.30]
    # cM0s = [1.0, 1.0, 1.0]

    assert(len(fA0s) == len(cM0s) and len(fA0s) == num_conditions)
    
    condition_ids = [f'c_{i}' for i in range(num_conditions)]
    conditions = [(cI0[i], fA0s[i]*cM0s[i], (1-fA0s[i])*cM0s[i]) 
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
    conditions_df = pd.DataFrame(columns=[CONDITION_ID, CONDITION_NAME, 'I', 'A', 'B'])
    
    for i, (cI0, cA0, cB0) in enumerate(conditions):
        # print(i, cI0, cA0, cB0)
        
        meas_a, meas_b = run_amici_simulation(amici_model, timepoints, cI0=cI0, cA0=cA0, cB0=cB0, sigma=meas_sigma)
        conv_a = (cA0 - meas_a)/cA0
        conv_b = (cB0 - meas_b)/cB0
        
        fA0 = cA0 / (cA0 + cB0)
        conv = conv_a * fA0 + conv_b * (1 - fA0)
        
        df = pd.DataFrame({
            # OBSERVABLE_ID
            SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
            TIME: timepoints,
            # MEASUREMENT: conv,
            'xA': conv_a,
            'xB': conv_b,
            'x': conv
        })
        
        if real:
            df_a = pd.DataFrame({
                OBSERVABLE_ID: ['obs_xa'] * len(timepoints),
                SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
                TIME: timepoints,
                MEASUREMENT: conv_a
            })
            
            df_b = pd.DataFrame({
                OBSERVABLE_ID: ['obs_xb'] * len(timepoints),
                SIMULATION_CONDITION_ID: [condition_ids[i]] * len(timepoints),
                TIME: timepoints,
                MEASUREMENT: conv_b
            })

        conditions_data = pd.DataFrame({
            CONDITION_ID: [condition_ids[i]],
            CONDITION_NAME: [condition_ids[i]],
            'I': [cI0],
            'A': [cA0],
            'B': [cB0]
        })
        
        if i == 0:
            if real:
                measurement_df = pd.concat([df_a, df_b], ignore_index=True)
            else:
                measurement_df = pd.concat([df], ignore_index=True)
            conditions_df = pd.concat([conditions_data], ignore_index=True)
        else:
            if real:
                measurement_df = pd.concat([measurement_df, df_a, df_b], ignore_index=True)
            else:
                measurement_df = pd.concat([measurement_df, df], ignore_index=True)
            conditions_df = pd.concat([conditions_df, conditions_data], ignore_index=True)
        
        
    return conditions_df, measurement_df

# problem --------------------------------------------------------------------

def get_petab_inputs(amici_model, sbml_model_filepath):
    
    parameters_df = define_parameters()
    observables_df, conditions_df, measurements_df = define_FRP_measurements(amici_model)
    
    return parameters_df, observables_df, conditions_df, measurements_df

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