from calendar import c
import numpy as np
import pandas as pd
# import petab
from petab.v1.C import PARAMETER_ID, PARAMETER_SCALE, LOWER_BOUND, UPPER_BOUND, NOMINAL_VALUE, ESTIMATE, OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT, OBSERVABLE_FORMULA, NOISE_FORMULA, LOG, CONDITION_ID, CONDITION_NAME, FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES
# import amici
from libsbml import Model
# import pypesto

# import PyPESTO.FRP.sbml as sbml
# import os
from PyPESTO.FRP.petab_ import plot_measurements

# from Sensitivity.FRP_Model_ import CopolymerizationModel
from PyPESTO.FRP import create_FRP2_v3
MODEL_NAME = 'FRP2_v3'

# Reactivity ratios
create_FRP2_v3.rA_true = rA_true = 0.5   # kpAA_true / kpAB_true
create_FRP2_v3.rB_true = rB_true = 18    # kpBB_true / kpBA_true
create_FRP2_v3.rX_true = rX_true = 1.0   # kpAA_true / kpBB_true?
create_FRP2_v3.KBB_true = KBB_true = 0.5  # kdBB_true / kpBB_true

# Initiator dissociation rate constant
create_FRP2_v3.kd_true = kd_true = 1e-4

# Propagation rate constant
create_FRP2_v3.kpAA_true = kpAA_true = 1e3
# kpAA_true = 1e3

# Termination rate constant
create_FRP2_v3.kt_true = kt_true = 1e4

kp_kt_ratio_true = kpAA_true / kt_true
kd_kt_true = kd_true * kt_true

# sbml_model_filepath = '/SBML/PyPESTO/FRP/sbml_model.xml'
amici_model, sbml_model_filepath = create_FRP2_v3.load_amici_from_sbml()
yaml_filepath = create_FRP2_v3.write_petab_files(amici_model, sbml_model_filepath)
observables_df, conditions_df, measurements_df = create_FRP2_v3.define_FRP_measurements(amici_model)



fig, axs = plot_measurements(measurements_df)
fig.savefig('FRP2_v3_measurements.png')

