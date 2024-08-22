import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import PyPESTO.FRP.sbml as sbml


# Define PEtab paths
FRP_model_dir = '/SBML/PyPESTO/FRP/FRP_model/'
FRP_model_xml_path = os.path.join(FRP_model_dir, 'FRP_model.xml')

FRP_model_conditions_path = os.path.join(FRP_model_dir, 'conditions.tsv')
FRP_model_measurements_path = os.path.join(FRP_model_dir, 'measurements.tsv')

# Create FRP model
model, document = sbml.create_FRP_model_v1()
sbml.outputSBML(document, FRP_model_xml_path)


timepoints = np.linspace(0, 100, 10)
fA0s = [0.25, 0.5, 0.75]
cM0s = [3.0, 3.0, 2.0]
M0 = 3.0
conditions = [(fA0*cM0[i], (1-fA0)*cM0[i]) for i, fA0 in enumerate(fA0s)]
num_conditions = len(conditions)
times = np.linspace(0, 3, 10)
sigma = 0.00

# Initialize an empty dataframe
measurement_df = pd.DataFrame(columns=['observableId', 'simulationConditionId', 'time', 'measurement'])

conditions_data = []

for i, (cA0, cB0) in enumerate(conditions, start=1):
    
    model.setInitialStates([5e-3, 0, cA0, cB0, 0, 0, 0, 0, 0, 0])
    
    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)
    
    meas_a = rdata.by_id('A') + sigma * np.random.randn(len(timepoints))
    meas_b = rdata.by_id('B') + sigma * np.random.randn(len(timepoints))
    
    # Create dataframes for obs_a and obs_b
    df_a = pd.DataFrame({
        'observableId': ['obs_a'] * len(timepoints),
        'simulationConditionId': [f'c_{i}'] * len(timepoints),
        'time': timepoints,
        'measurement': meas_a
    })

    df_b = pd.DataFrame({
        'observableId': ['obs_b'] * len(timepoints),
        'simulationConditionId': [f'c_{i}'] * len(timepoints),
        'time': timepoints,
        'measurement': meas_b
    })

    # Append to the main dataframe
    measurement_df = pd.concat([measurement_df, df_a, df_b], ignore_index=True)
    
    # Add to conditions data
    conditions_data.append({
                            'conditionId': f'c_{i}', 
                            'conditionName': f'c_{i}',
                            'A': cA0,
                            'B': cB0})

# Convert conditions data to dataframe
conditions_df = pd.DataFrame(conditions_data)

print(conditions_df)
print(measurement_df)
# Save the dataframe to a file
measurement_df.to_csv('/SBML/PyPESTO/FRP/FRP_model/multiple_conditions/measurements.tsv', sep='\t', index=False)
conditions_df.to_csv('/SBML/PyPESTO/FRP/FRP_model/multiple_conditions/conditions.tsv', sep='\t', index=False)
# measurement_df.to_csv('doc/example/TEST_conversion_reaction/multiple_conditions/measurements.tsv', sep='\t', index=False)
# measurement_df.to_csv('/SBML/PyPESTO/TEST_conversion/TEST_conversion_reaction/measurements.tsv', sep='\t', index=False)

# conditions_df.to_csv('/SBML/PyPESTO/TEST_conversion/TEST_conversion_reaction/multiple_conditions/conditions.tsv', sep='\t', index=False)