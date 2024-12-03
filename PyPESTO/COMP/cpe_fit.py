from PyPESTO.COMP.CPE import CopolymerEquations, cpe_model, create_wrapper
import matplotlib.pyplot as plt
import numpy as np

# Define the model parameters
r1_true = 3.0
r2_true = 0.5
b1_true = 0.0
g1_true = 0.0
b2_true = 0.0
g2_true = 0.0

# Initial fit guesses
r1_guess = r1_true
r2_guess = r2_true

params_true = {
    'r1': r1_true,
    'r2': r2_true,
    'b1': b1_true,
    'g1': g1_true,
    'b2': b2_true,
    'g2': g2_true,
}

params_fit_guess = {
    'r1': r1_guess,
    'r2': r2_guess,
}

params_fit_bounds = {
    'r1': (1e-3, 1e3),
    'r2': (1e-3, 1e3),
}

# Assert that params_fit is within params_fit_bounds
if not all([params_fit_guess[key] >= params_fit_bounds[key][0] and params_fit_guess[key] <= params_fit_bounds[key][1] for key in params_fit_guess]):
    raise ValueError('Initial guess is not within bounds')

params_fixed = {key:val for key, val in params_true.items() if key not in params_fit_guess}


print('True param values:', params_true)
print('Fixed params:', params_fixed)
print('Fit params (bounds):', params_fit_bounds)
print('Fit params (initial guess):', params_fit_guess)

X_eval = np.linspace(0.02, 0.95, 10)
data = cpe_model(X_eval, **params_true)

rng = np.random.default_rng()
y_noise = 0.00 * rng.normal(size=data.size)
y_exp = data + y_noise

# Plot experimental data
data_exp = y_exp.reshape(2, -1)
plt.figure(dpi=150)
# plt.plot(X_eval, data[0], 'o-')
# plt.plot(X_eval, 2*X_eval - data[1], 'o-')
plt.plot(data_exp[0], data_exp[1], 'o-')
plt.plot([0, 1], [0, 1], 'k-', alpha=0.15)
plt.xlim([0, 1])
plt.ylim([0, 1])

# Make the plot square
plt.gca().set_aspect('equal', adjustable='box')

# Increase the font size
plt.rcParams.update({'font.size': 16})

plt.xlabel('A conversion')
plt.ylabel('B conversion')

# Parameter estimation
from scipy.optimize import curve_fit
fit_func = create_wrapper(cpe_model, params_fixed, params_fit_guess)
params, cov = curve_fit(fit_func, X_eval, y_exp, method='trf', loss='soft_l1', p0=list(params_fit_guess.values()), xtol=1e-15, ftol=1e-15)#, bounds=params_fit_bounds.values(), verbose=2, full_output=True, xtol=1e-15, ftol=1e-15, gtol=1e-15)
