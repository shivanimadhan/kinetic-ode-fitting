
import numpy as np
import os
import amici

MODEL_NAME = 'FRP2_v4'

def cpe_sim(model: amici.AmiciModel, X_eval, fA0, log_r=None):

    # Set parameter values
    if log_r is None:
        rA = 18.0
        rB = 0.5
        rX = 1.0
        KAA = 0.5
        KAB = 0.0
        KBA = 0.0
        KBB = 0.0
    else:
        log_rA, log_rB, log_rX, KAA, KAB, KBA, KBB = log_r
        rA = 10.**log_rA
        rB = 10.**log_rB
        rX = 10.**log_rX
        
    model.setParameterByName('kpAA', 1)
    model.setParameterByName('rA', rA)
    model.setParameterByName('rB', rB)
    model.setParameterByName('rX', rX)
    model.setParameterByName('KAA', KAA)
    model.setParameterByName('KAB', KAB)
    model.setParameterByName('KBA', KBA)
    model.setParameterByName('KBB', KBB)
    
    # Ensure sensitivities are computed for parameters in log_r
    # parameter_ids = ['log_rA', 'log_rB', 'log_rX', 'KAA', 'KAB', 'KBA', 'KBB']
    # model.setParameterList(parameter_ids)
    
    parameter_ids = model.getParameterIds()
    param_indices = {param: idx for idx, param in enumerate(parameter_ids)}
    
    model.requireSensitivitiesForAllParameters()
    
    # Solver settings
    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-12)
    solver.setRelativeTolerance(1e-8)
    solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    # solver.setSensitivityRelativeTolerance(1e-8)
    # solver.setSensitivityAbsoluteTolerance(1e-12)
    
    # Set initial conditions and timepoints
    model.setInitialStates([0.0, 0.005, fA0, 1-fA0, 0, 0, 0, 0, 0, 0, 0])
    model.setTimepoints(np.linspace(0, 4000, 4000))
    
    # Run simulation
    rdata = amici.runAmiciSimulation(model, solver)
    
    return rdata  # rdata contains both outputs and sensitivities

def get_cpe_output(rdata, X_eval, model):
    import numpy as np

    # Get observable IDs
    observable_ids = model.getObservableIds()
    # idx_A = observable_ids.index('A')
    # idx_B = observable_ids.index('B')

    # Extract outputs
    cA = rdata['y'][:, 0]
    cB = rdata['y'][:, 1]
    
    # print('cA', cA.shape)

    # Extract sensitivities
    sens_cA = rdata['sy'][:, :, 0]  # Shape: (4000, 10)
    sens_cB = rdata['sy'][:, :, 1]  # Shape: (4000, 10)

    # Compute xA and xB
    cA0 = cA[0]
    cB0 = cB[0]
    xA = 1 - cA / cA0
    xB = 1 - cB / cB0
    x = 1 - (cA + cB) / (cA0 + cB0)

    # Compute sensitivities of xA and xB
    sens_xA = -sens_cA / cA0
    sens_xB = -sens_cB / cB0
    
    # print('sens_xA', sens_xA.shape)

    # Adjust sensitivities using chain rule
    parameter_ids = model.getParameterIds()
    param_indices = {param: idx for idx, param in enumerate(parameter_ids)}

    # Create adjustments for parameters
    adjustments = np.ones(len(parameter_ids))
    for param in ['rA', 'rB', 'rX']:
        idx = param_indices[param]
        param_value = model.getParameterByName(param)
        adjustments[idx] = np.log(10) * param_value

    # print('adjustments', adjustments.shape)
    # Apply adjustments to sensitivities
    sens_xA *= adjustments
    sens_xB *= adjustments

    # Select sensitivities for estimated parameters
    estimated_params = ['rA', 'rB', 'rX', 'KAA', 'KAB', 'KBA', 'KBB']
    estimated_indices = [param_indices[param] for param in estimated_params]
    sens_xA = sens_xA[:, estimated_indices]
    sens_xB = sens_xB[:, estimated_indices]

    # Handle interpolation for outputs
    xA_interp = np.interp(X_eval, x, xA)
    xB_interp = np.interp(X_eval, x, xB)

    # Interpolate sensitivities
    sens_xA_interp = np.zeros((len(X_eval), len(estimated_indices)))
    sens_xB_interp = np.zeros((len(X_eval), len(estimated_indices)))
    for i in range(len(estimated_indices)):
        sens_xA_interp[:, i] = np.interp(X_eval, x, sens_xA[:, i])
        sens_xB_interp[:, i] = np.interp(X_eval, x, sens_xB[:, i])

    return xA_interp, xB_interp, sens_xA_interp, sens_xB_interp


# def get_cpe_output(rdata: amici.ReturnData, X_eval):
#     # cA = rdata['y'][:, rdata['y_names'].index('A')]
#     # cB = rdata['y'][:, rdata['y_names'].index('B')]
    
    
#     cA = rdata.by_id('A')
#     cB = rdata.by_id('B')
    
#     # Sensitivities
#     sens_cA = rdata['sy'][:, 0, :]
#     sens_cB = rdata['sy'][:, 1, :]
#     # sens_cA = rdata['sy'][:, rdata['y_names'].index('A'), :]
#     # sens_cB = rdata['sy'][:, rdata['y_names'].index('B'), :]
    
#     # Compute xA and xB
#     cA0 = cA[0]
#     cB0 = cB[0]
#     xA = 1 - cA / cA0
#     xB = 1 - cB / cB0
#     x  = 1 - (cA + cB) / (cA0 + cB0)
    
#     # Compute sensitivities of xA and xB
#     # Using chain rule
#     sens_xA = - sens_cA / cA0
#     sens_xB = - sens_cB / cB0
    
#     # Handle interpolation
#     x_interp = x
#     xA_interp = xA
#     xB_interp = xB
#     # Interpolated outputs
#     xA = np.interp(X_eval, x_interp, xA_interp)
#     xB = np.interp(X_eval, x_interp, xB_interp)
    
#     # Interpolate sensitivities
#     sens_xA_interp = np.empty((len(X_eval), sens_xA.shape[1]))
#     sens_xB_interp = np.empty((len(X_eval), sens_xB.shape[1]))
#     for i in range(sens_xA.shape[1]):  # Loop over parameters
#         sens_xA_interp[:, i] = np.interp(X_eval, x_interp, sens_xA[:, i])
#         sens_xB_interp[:, i] = np.interp(X_eval, x_interp, sens_xB[:, i])
    
#     return xA, xB, sens_xA_interp, sens_xB_interp


def cpe_model(X_eval, log_r=None):

    model_output_dir = os.path.join('tmp/', MODEL_NAME)
    model_module = amici.import_model_module(MODEL_NAME, model_output_dir)
    model = model_module.getModel()
    
    # List to collect outputs and sensitivities
    fit_data_list = []
    sens_fit_data_list = []
    
    initial_fA0s = [0.25, 0.5, 0.75]
    for fA0 in initial_fA0s:
        rdata = cpe_sim(model, X_eval, fA0, log_r)
        xA, xB, sens_xA, sens_xB = get_cpe_output(rdata, X_eval, model)
        
        # Combine xA and xB
        fit_data = np.concatenate((xA, xB))
        fit_data_list.append(fit_data)
        
        # Combine sensitivities
        sens_fit_data = np.vstack((sens_xA, sens_xB))
        sens_fit_data_list.append(sens_fit_data)
    
    # Combine all data
    fit_data = np.concatenate(fit_data_list)
    sens_fit_data = np.vstack(sens_fit_data_list)
    
    return fit_data, sens_fit_data  # Return both outputs and sensitivities


# def acqf_FRP2_v4(log_r):

#     X_eval = np.linspace(0.1, 0.8, 8)
#     exp_data = np.array([0.15853124, 0.30176469, 0.43230576, 0.55115553, 0.65896751,
#        0.75597703, 0.84182983, 0.91513637, 0.08048959, 0.16607844,
#        0.25589808, 0.34961482, 0.44701083, 0.54800766, 0.65272339,
#        0.76162121, 0.14445154, 0.27676588, 0.39640375, 0.50484807,
#        0.6037441 , 0.69456924, 0.77856031, 0.85674143, 0.05554846,
#        0.12323412, 0.20359625, 0.29515193, 0.3962559 , 0.50543076,
#        0.62143969, 0.74325857, 0.12385877, 0.24530836, 0.36031259,
#        0.46568158, 0.56064598, 0.64670672, 0.72591171, 0.77765358,
#        0.02842368, 0.06407492, 0.11906224, 0.20295527, 0.31806207,
#        0.45987985, 0.62226486, 0.74309089])
#     fit_data = cpe_model(X_eval, log_r=log_r)
    
#     ssr = np.sum((exp_data - fit_data)**2)
#     return np.log(ssr)

def acqf_FRP2_v4_wrapper(exp_data):

    def acqf_FRP2_v4(log_r):
        
        X_eval = np.linspace(0.1, 0.8, 8)
        # exp_data = np.array([0.15853124, 0.30176469, 0.43230576, 0.55115553, 0.65896751,
        #    0.75597703, 0.84182983, 0.91513637, 0.08048959, 0.16607844,
        #    0.25589808, 0.34961482, 0.44701083, 0.54800766, 0.65272339,
        #    0.76162121, 0.14445154, 0.27676588, 0.39640375, 0.50484807,
        #    0.6037441 , 0.69456924, 0.77856031, 0.85674143, 0.05554846,
        #    0.12323412, 0.20359625, 0.29515193, 0.3962559 , 0.50543076,
        #    0.62143969, 0.74325857, 0.12385877, 0.24530836, 0.36031259,
        #    0.46568158, 0.56064598, 0.64670672, 0.72591171, 0.77765358,
        #    0.02842368, 0.06407492, 0.11906224, 0.20295527, 0.31806207,
        #    0.45987985, 0.62226486, 0.74309089])
        
        # exp_data = np.array([0.14401001, 0.3093343 , 0.4305011 , 0.54237473, 0.6461774 ,
        # 0.75351303, 0.85821263, 0.92227847, 0.07456658, 0.16950477,
        # 0.25025108, 0.34987644, 0.44769104, 0.5520851 , 0.64218282,
        # 0.77815377, 0.1383821 , 0.27454378, 0.39647061, 0.51328827,
        # 0.60196285, 0.69827501, 0.76301668, 0.85529159, 0.06181636,
        # 0.12305216, 0.22687549, 0.29945819, 0.39514776, 0.50757883,
        # 0.63516968, 0.73359255, 0.12389453, 0.2447328 , 0.34686853,
        # 0.45213986, 0.54638493, 0.64316831, 0.72316651, 0.77838068,
        # 0.02430458, 0.05946684, 0.12011565, 0.20031807, 0.32114777,
        # 0.46153402, 0.62394588, 0.74379217])
        
        # Ensure log_r is a NumPy array
        log_r = np.array(log_r)
        
        # Obtain simulated data and sensitivities
        fit_data, sens_fit_data = cpe_model(X_eval, log_r=log_r)
        
        # Compute residuals
        residuals = exp_data - fit_data
        
        # Compute SSR
        ssr = np.sum(residuals**2)
        
        # Handle potential division by zero
        if ssr == 0:
            ssr = 1e-12  # Small value to prevent division by zero
        
        # Compute gradient
        # Shape of residuals: (N,)
        # Shape of sens_fit_data: (N, num_parameters)
        gradient = - (2 / ssr) * np.dot(residuals, sens_fit_data)
        
        # Return the objective function value and gradient
        return (np.log(ssr), gradient)
    
    return acqf_FRP2_v4

    
    # exp_data = cpe_model(X_eval, log_r=None)
    
    # With noise
    # exp_data = np.array([0.18214712, 0.34695785, 0.49580307, 0.63839025, 0.76689604,
    #    0.82746461, 0.90594965, 0.07463086, 0.15339151, 0.21755202,
    #    0.31765773, 0.41612244, 0.51800079, 0.63912797, 0.13906993,
    #    0.28698877, 0.43560745, 0.56726073, 0.71564459, 0.79579776,
    #    0.89047681, 0.06685354, 0.09896311, 0.17456238, 0.23602791,
    #    0.32760088, 0.40004966, 0.51468394, 0.11918277, 0.24234121,
    #    0.35334949, 0.47404672, 0.59448669, 0.69999133, 0.81039383,
    #    0.05014606, 0.0838756 , 0.10909381, 0.16476532, 0.22325151,
    #    0.28177951, 0.3838005 ])
    
    # No noise
    # exp_data = np.array([0.18002763, 0.35072399, 0.50383562, 0.6378749 , 0.75153614,
    #    0.84384354, 0.91432327, 0.07332412, 0.14975867, 0.23205479,
    #    0.32070837, 0.41615462, 0.51871882, 0.62855891, 0.14851707,
    #    0.29506282, 0.43504387, 0.56672048, 0.68770204, 0.79473025,
    #    0.88355261, 0.05148293, 0.10493718, 0.16495613, 0.23327952,
    #    0.31229796, 0.40526975, 0.51644739, 0.12067434, 0.24138814,
    #    0.36060751, 0.4778742 , 0.59247597, 0.70322011, 0.80789916,
    #    0.03797697, 0.07583557, 0.11817747, 0.16637741, 0.22257208,
    #    0.29033966, 0.37630252])
    
    