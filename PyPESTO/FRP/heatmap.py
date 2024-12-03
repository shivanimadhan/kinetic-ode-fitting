# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

# from PyPESTO.FRP import FRP2_v3

# # x axis: kp_kt_ratio
# # y axis: kd_kt

# # Create a heatmap of some kind of measure of data "quality"
# # Bad data:
# # - data stays at initial value (too slow)
# # - data instantly goes to zero (too fast)

# # Create a grid of parameter values
# kp_kt_ratio_values = np.logspace(-6, 6, 100, base=10)
# kd_kt_values = np.logspace(-6, 6, 100, base=10)

# # Run the simulations for these parameters
# # Calculate the data quality metric for each set of parameters (sum of variances)
# # Plot the heatmap

# # kd, kp, kt

from calendar import c
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter, LogLocator

from petab.v1.C import PARAMETER_ID, PARAMETER_SCALE, LOWER_BOUND, UPPER_BOUND, NOMINAL_VALUE, ESTIMATE, OBSERVABLE_ID, SIMULATION_CONDITION_ID, TIME, MEASUREMENT, OBSERVABLE_FORMULA, NOISE_FORMULA, LOG, CONDITION_ID, CONDITION_NAME, FORMAT_VERSION, PARAMETER_FILE, PROBLEMS, SBML_FILES, CONDITION_FILES, MEASUREMENT_FILES, OBSERVABLE_FILES
from libsbml import Model

from PyPESTO.FRP.petab_ import plot_measurements

from PyPESTO.FRP import create_heatmap
MODEL_NAME = 'HEATMAP'

def sum_of_variances(measurements_df):
    conv_df = measurements_df.copy()
    conditions = measurements_df[SIMULATION_CONDITION_ID].unique()
    observables = measurements_df[OBSERVABLE_ID].unique()
    variance_dict = {}

    for c in conditions:
        for o in observables:
            meas_data = conv_df.loc[(conv_df[SIMULATION_CONDITION_ID] == c) & (conv_df[OBSERVABLE_ID] == o), MEASUREMENT]
            max_conc = meas_data.max()
            meas_data = conv_df.loc[(conv_df[SIMULATION_CONDITION_ID] == c) & (conv_df[OBSERVABLE_ID] == o) & (conv_df[TIME] > 0), MEASUREMENT]
            conv_data = (max_conc - meas_data) / max_conc
        
            variance_dict[(c, o)] = conv_data.var()

    return np.sum(list(variance_dict.values()))

def generate_heatmap(variances, kp_kt_ratio_values, kd_kt_values):
    df = pd.DataFrame(variances, index=kp_kt_ratio_values, columns=kd_kt_values)

    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(df, annot=False, cmap='viridis', cbar=True)

    plt.xlabel('kp_kt_ratio')
    plt.ylabel('kd_kt')
    plt.title('Heatmap of "Quality" of Data')
    
    file_path = f'/SBML/PyPESTO/FRP/{MODEL_NAME}/heatmap.png'
    plt.savefig(file_path, bbox_inches='tight')

    plt.close()

# Reactivity ratios
create_heatmap.rA_true = rA_true = 0.5   # kpAA_true / kpAB_true
create_heatmap.rB_true = rB_true = 18    # kpBB_true / kpBA_true
create_heatmap.rX_true = rX_true = 1.0   # kpAA_true / kpBB_true?
create_heatmap.KBB_true = KBB_true = 0.5  # kdBB_true / kpBB_true

# Initiator dissociation rate constant
create_heatmap.kd_true = kd_true = 1e-4

# Propagation rate constant
create_heatmap.kpAA_true = kpAA_true = 1e3

# Termination rate constant
create_heatmap.kt_true = kt_true = 1e4

# Sampling values for kp_kt_ratio, kd_kt, and kpAA
num_points = 100
kp_kt_ratio_values = np.logspace(-6, 6, num_points, base=10)
kd_kt_values = np.logspace(-6, 6, num_points, base=10)
#kpAA_values = np.logspace(-6, 6, 100, base=10) 

variances = np.zeros((num_points, num_points))

amici_model, sbml_model_filepath = create_heatmap.load_amici_from_sbml()

for i, point_1 in enumerate(kp_kt_ratio_values):
    create_heatmap.kp_kt_ratio_true = point_1

    for j, point_2 in enumerate(kd_kt_values):
        create_heatmap.kd_kt_true = point_2

        observables_df, conditions_df, measurements_df = create_heatmap.define_FRP_measurements(amici_model)
        
        # Debugging prints
        # print(f"kp_kt_ratio: {point_1}, kd_kt: {point_2}")
        # print(f"measurements_df.head():\n{measurements_df.head()}")

        metric = sum_of_variances(measurements_df)
        variances[i][j] = metric
        # print(f"Metric: {metric}")


generate_heatmap(variances, kp_kt_ratio_values, kd_kt_values)