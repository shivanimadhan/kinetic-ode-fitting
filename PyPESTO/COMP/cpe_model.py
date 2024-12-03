
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import plotly.graph_objects as go
import amici

MODEL_NAME = 'FRP2_v4'

def cpe_sim(model, X_eval, fA0, log_r=None):
    
    # rA, rB, rX, KAA, KAB, KBA, KBB, kpAA, kt_kp_ratio, kd_kt = 10**(log_r)
    if log_r is None:
        rA = 18.0
        rB = 0.5
        rX = 1.0
        KAA = 0.5
        KAB = 0.0
        KBA = 0.0
        KBB = 0.0
    # else:
    #     log_rA, log_rB, log_rX, KAA, KAB, KBA, KBB = log_r
    #     rA = 10.**log_rA
    #     rB = 10.**log_rB
    #     rX = 10.**log_rX
    else:
        log_rA, log_rB, log_rX, z_KAA, z_KAB, z_KBA, z_KBB = log_r
        rA = 10.**log_rA
        rB = 10.**log_rB
        rX = 10.**log_rX
        
        KAA = z_KAA**2
        KAB = z_KAB**2
        KBA = z_KBA**2
        KBB = z_KBB**2
        
    model.setParameterByName('kpAA', 1)
    model.setParameterByName('rA', rA)
    model.setParameterByName('rB', rB)
    model.setParameterByName('rX', rX)
    model.setParameterByName('KAA', KAA)
    model.setParameterByName('KAB', KAB)
    model.setParameterByName('KBA', KBA)
    model.setParameterByName('KBB', KBB)
    model.setParameterByName('kt', 0) #kt_kp_ratio * kpAA)
    model.setParameterByName('kd', 0) #kd_kt / kt_kp_ratio)
    # model.requireSensitivitiesForAllParameters()
    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-12)
    # solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    # solver.setSensitivityOrder(amici.SensitivityOrder_first)

    # model.setInitialStates([0.0, 0.005, fA0, 1-fA0, 0, 0, 0, 0, 0, 0])


    model.setInitialStates([0.0, 0.02, fA0, 1-fA0, 0, 0, 0, 0, 0, 0, 0])
    model.setTimepoints(np.linspace(0, 4000, 4000))
    rdata = amici.runAmiciSimulation(model, solver)
    
    return rdata

def get_cpe_output(rdata, X_eval):
    
    cA = rdata.by_id('A')
    cB = rdata.by_id('B')
    
    xA = 1 - cA / cA[0]
    xB = 1 - cB / cB[0]
    x  = 1 - (cA + cB) / (cA[0] + cB[0])
    
    xA = np.interp(X_eval, x, xA)
    xB = np.interp(X_eval, x, xB)
    
    return xA, xB

def cpe_model(X_eval, log_r=None, exp_noise=0.0, fA0s=[0.5]):

    model_output_dir = os.path.join('tmp/', MODEL_NAME)
    model_module = amici.import_model_module(MODEL_NAME, model_output_dir)
    model = model_module.getModel()
    
    fit_data = []
    
    for fA0 in fA0s:
        # Get outputs for the current fA0
        xA, xB = get_cpe_output(cpe_sim(model, X_eval, fA0, log_r), X_eval)
        
        # Add xA and xB to the fit_data in the required order
        fit_data.extend([xA, xB])

    # Combine into a single array while preserving the order
    fit_data = np.concatenate(fit_data)

    # Add noise if specified
    fit_data = fit_data + np.random.normal(0, exp_noise, fit_data.shape)
    
    return fit_data



def acqf_FRP2_v4_wrapper(exp_data, x_eval, fA0s):
    
    def acqf_FRP2_v4(log_r):
        fit_data = cpe_model(x_eval, log_r=log_r, fA0s=fA0s)

        ssr = np.sum((exp_data - fit_data)**2)
        return np.log(ssr)

    return acqf_FRP2_v4

import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_data(xdata, y_vals, ax=None, n_exp=3, plot_style="both", subplot=False):
    x = xdata.reshape(n_exp, 2, -1)
    
    # Define colormaps and specify a range within the colormap to avoid extreme colors
    cmap_red = cm.get_cmap('Reds')
    cmap_blue = cm.get_cmap('Blues')
    
    # Select colors within a specified range (e.g., 0.3 to 0.7) to avoid extremes
    colors_red = [cmap_red(0.3 + i * 0.2) for i in range(n_exp)]
    colors_blue = [cmap_blue(0.3 + i * 0.2) for i in range(n_exp)]
    
    # Define different marker styles for each dataset
    markers = ['o', 's', '^']  # Circle, square, and triangle markers

    if subplot:
        if ax is None:
            fig, axes = plt.subplots(1, n_exp, figsize=(5 * n_exp, 5), sharey=True)
            if n_exp > 1:
                axes = axes.flatten()
            else:
                axes = [axes]
        else:
            axes = ax
    elif ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
        axes = [ax]
    else:
        axes = [ax]
        
    if axes is None:
        return None

    # Plot each dataset
    for i in range(n_exp):
        current_ax = axes[i] if subplot else ax
        if plot_style == "both":
            current_ax.plot(y_vals, x[i][0], linestyle='-', color=colors_red[i], marker=markers[i], label=f"Data {i+1} Red" if i == 0 else "")
            current_ax.plot(y_vals, x[i][1], linestyle='-', color=colors_blue[i], marker=markers[i], label=f"Data {i+1} Blue" if i == 0 else "")
        elif plot_style == "markers":
            current_ax.scatter(y_vals, x[i][0], color=colors_red[i], marker=markers[i], label=f"Data {i+1} Red" if i == 0 else "")
            current_ax.scatter(y_vals, x[i][1], color=colors_blue[i], marker=markers[i], label=f"Data {i+1} Blue" if i == 0 else "")
        elif plot_style == "lines":
            current_ax.plot(y_vals, x[i][0], linestyle='-', color=colors_red[i], label=f"Data {i+1} Red" if i == 0 else "")
            current_ax.plot(y_vals, x[i][1], linestyle='-', color=colors_blue[i], label=f"Data {i+1} Blue" if i == 0 else "")

        current_ax.plot([0, 1], [0, 1], 'k--', alpha=0.15)

        current_ax.set_xlim(0, 1)
        current_ax.set_ylim(0, 1)
        current_ax.set_xlabel('Total Conversion')
        # if i == 0 or not subplot:
        current_ax.set_ylabel('Monomer Conversion')
        current_ax.legend(['A', 'B'], loc='upper left', frameon=False)  # Add legend if desired
    
    
    # if subplot:
    #     axes[0].legend()  # Add legend to the first subplot if desired
    # else:
    #     ax.legend()  # Add legend if desired

    return axes if subplot else ax



def transform_obj_for_plot(obj_df):
    
    # Pivot the DataFrame to create matrices for plotting
    z_matrix = obj_df.pivot(index='y', columns='x', values='z').values
    x_unique = np.sort(obj_df['x'].unique())
    y_unique = np.sort(obj_df['y'].unique())

    return z_matrix, x_unique, y_unique

def plot_acqf_surface(obj_df, xaxis, yaxis, **kwargs):
    # Pivot the DataFrame to create matrices for plotting
    z_matrix, x_unique, y_unique = transform_obj_for_plot(obj_df)

    # Create the 3D surface plot with hover text
    surface = go.Figure(data=[go.Surface(
        z=z_matrix,
        x=x_unique,
        y=y_unique,
        hovertemplate=f"<b>{xaxis}</b>: %{{x}}<br><b>{yaxis}</b>: %{{y}}<br><b>Acquisition Function</b>: %{{z}}<extra></extra>"
    )])

    surface.update_layout(
        scene=dict(
            xaxis_title=xaxis,
            yaxis_title=yaxis,
            zaxis_title='Acquisition Function',
        ),
        width=800,
        height=600,
        **kwargs
    )
    surface.show()

def plot_acqf_contour(obj_df, xaxis, yaxis, **kwargs):
    # Pivot the DataFrame to create matrices for plotting
    z_matrix, x_unique, y_unique = transform_obj_for_plot(obj_df)

    # Create the 2D contour plot with hover text
    contour = go.Figure(data=[go.Contour(
        z=z_matrix,
        x=x_unique,
        y=y_unique,
        hovertemplate=f"<b>{xaxis}</b>: %{{x}}<br><b>{yaxis}</b>: %{{y}}<br><b>Acquisition Function</b>: %{{z}}<extra></extra>"
    )])

    contour.update_layout(
        xaxis_title=xaxis,
        yaxis_title=yaxis,
        width=800,
        height=600,
        **kwargs
    )
    contour.show()
