import numpy as np
from scipy.optimize import fsolve, least_squares
from scipy.integrate import odeint, solve_ivp
import amici

class CopolymerEquations:
    
    def __init__(self, inputs):
        
        self.fA0 = inputs[0]
        self.fB0 = 1 - self.fA0
        
        self.kpAA = self.k1 = inputs[1]
        self.kdAA = self.k2 = inputs[2]
        self.kpAB = self.k3 = inputs[3]
        self.kdAB = self.k4 = inputs[4]
        self.kpBA = self.k5 = inputs[5]
        self.kdBA = self.k6 = inputs[6]
        self.kpBB = self.k7 = inputs[7]
        self.kdBB = self.k8 = inputs[8]
        
        self.bounds = [0, 0.9]
        self.num_points = 100
        self.t_eval = np.linspace(self.bounds[0], self.bounds[1], self.num_points)
        self.method = 'BDF'
        
    def set_bounds(self, bounds):
        self.bounds = bounds
        self.t_eval = np.linspace(self.bounds[0], self.bounds[1], self.num_points)
        
    @staticmethod
    def fromLundberg(inputs):
        
        f1, r1, r2, b1, g1, b2, g2 = inputs
        new_inputs = [
            f1,
            1,
            b1,
            1/r1,
            g2,
            1/r2,
            g1,
            1,
            b2
        ]
        
        return CopolymerEquations(new_inputs)

    def solve_izu(self, X_eval=None):
        
        def izu_vars(X, A, B):
            def izu_var_eqns(vars):
                a, b, c, d = vars
                eq1 = a + b - 1
                eq2 = a*(self.k3*B + (1-c)*self.k6) - b*(self.k5*A+(1-d)*self.k4)
                eq3 = b*c*(1-d)*self.k4 - a*(c*(self.k1*A + self.k2 + self.k3*B) - (self.k1*A + c*c*self.k2))
                eq4 = a*d*(1-c)*self.k6 - b*(d*(self.k7*B + self.k8 + self.k5*A) - (self.k7*B + d*d*self.k8))
                return [eq1, eq2, eq3, eq4]
            initial_guess = [0.5, 0.5, 0.5, 0.5]
            res = least_squares(izu_var_eqns, initial_guess, bounds=([0, 0, 0, 0], [1, 1, 1, 1]))
            # print(X, res.cost, res.x)
            return res.x
            # initial_guess = [0.5, 0.5, 0.5, 0.5]
            # return fsolve(witmer_var_eqns, initial_guess,xtol=1e-10, maxfev=1000)
            
        def izu_eqns(X, fA):
            fA = fA[0]
            A = fA*(1-X)*(self.fA0 + self.fB0)
            B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
            a,b,c,d = izu_vars(X, A, B)

            dAdt = a*self.k1*A + b*self.k5*A - a*((1-c)*self.k6 + c*self.k2)
            dBdt = b*self.k7*B + a*self.k3*B - b*((1-d)*self.k4 + d*self.k8)
            FA = dAdt / (dAdt + dBdt)
            
            return [(fA - FA)/(1-X)]
        
        sol = solve_ivp(izu_eqns, self.bounds, [self.fA0], method=self.method, max_step=0.01, t_eval=self.t_eval)
        # return sol
        # X = sol.t
        # fA = sol.y[0]
        
        # A = fA*(1-X)*(self.fA0 + self.fB0)
        # B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
        
        if X_eval is None:
            X_eval = sol.t
            fA = sol.y[0]
        else:
            # interpolate fA for X_eval
            fA = np.interp(X_eval, sol.t, sol.y[0])
        #     A = fA*(1-X_eval)*(self.fA0 + self.fB0)
        #     B = (1 - fA)*(1-X_eval)*(self.fA0 + self.fB0)
        # else:
        #     A = fA*(1-X)*(self.fA0 + self.fB0)
        #     B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
        
        A = fA*(1-X_eval)*(self.fA0 + self.fB0)
        B = (1 - fA)*(1-X_eval)*(self.fA0 + self.fB0)
        
        xA = (self.fA0 - A)/(self.fA0)
        xB = (self.fB0 - B)/(self.fB0)
        
        return X_eval, xA, xB
    
    def solve_wittmer(self):
        
        q1 = self.kdAB/self.kpBA
        q2 = self.kdBA/self.kpAB
        r1 = self.kpAA/self.kpAB 
        r2 = self.kpBB/self.kpBA
        K1 = self.kdAA/self.kpAA
        K2 = self.kdBB/self.kpBB
        
        def wittmer_vars(X, A, B):
            def wittmer_var_eqns(vars):
                x1, y1 = vars
                eq1 = q1*y1*(B+q2*x1)/(A+q1*y1) - B - r1*(A+K1)+2*r1*(A+K1*(1-x1)**2)/(2*(1-x1))
                eq2 = q2*x1*(A+q1*y1)/(B+q2*x1) - A - r2*(B+K2)+2*r2*(B+K2*(1-y1)**2)/(2*(1-y1))
                return [eq1, eq2]
            initial_guess = [0.5, 0.5]
            # return fsolve(wittmer_var_eqns, initial_guess,xtol=1e-10, maxfev=1000)
            res = least_squares(wittmer_var_eqns, initial_guess, bounds=([0, 0], [1, 1]))
            # print(X, res.cost, res.x)
            return res.x
            
        def wittmer_eqns(X, fA):
            fA = fA[0] + 1e-40
            X += 1e-40
            A = fA*(1-X)*(self.fA0 + self.fB0)
            B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
            x1, y1 = wittmer_vars(X, A, B)
            
            dAdt = 1 + (r1*A/B - r1*K1*(1-x1)/B) / (1 - q1*y1/B*(B+q2*x1)/(A+q1*y1))
            dBdt = 1 + (r2*B/A - r2*K2*(1-y1)/A) / (1 - q2*x1/A*(A+q1*y1)/(B+q2*x1))
            FA = dAdt / (dAdt + dBdt)
            return [(fA - FA)/(1-X)]
        
        sol = solve_ivp(wittmer_eqns, self.bounds, [self.fA0], method=self.method, max_step=0.001, t_eval=self.t_eval)
        X = sol.t
        fA = sol.y[0]
        
        A = fA*(1-X)*(self.fA0 + self.fB0)
        B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
        
        xA = (self.fA0 - A)/(self.fA0)
        xB = (self.fB0 - B)/(self.fB0)
        
        return X, xA, xB
        
    def solve_kruger(self):
        
        rA  = self.k1 / self.k3 
        rB  = self.k7 / self.k5
        RA  = self.k4 / self.k5
        RAA = self.k2 / self.k3
        RB  = self.k6 / self.k3
        RBB = self.k8 / self.k5
        
        def kruger_vars(X, A, B):
            def kruger_var_eqns(vars):
                PAA, PAB, PBA, PBB = vars
                a = 1 - RA*PBA / (A + RA*PBA)
                b = RAA - RA*RB*PBA / (A + RA*PBA)
                c = 1 - RB*PAB / (B + RB*PAB)
                d = RBB - RB*RA*PAB / (B + RB*PAB)
                
                eq1 = b*PAB*PAB + (rA*A + B*a - b)*PAB - B*a
                eq2 = d*PBA*PBA + (rB*B + A*c - d)*PBA - A*c
                eq3 = PAA + PAB - 1
                eq4 = PBB + PBA - 1
                
                return [eq1, eq2, eq3, eq4]
            initial_guess = [0.5, 0.5, 0.5, 0.5]
            res = least_squares(kruger_var_eqns, initial_guess, bounds=([0, 0, 0, 0], [1, 1, 1, 1]))
            # print(X, res.cost, res.x)
            return res.x
            
        def kruger_eqns(X, fA):
            fA = fA[0] + 1e-40
            X += 1e-40
            A = fA*(1-X)*(self.fA0 + self.fB0)
            B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
            PAA, PAB, PBA, PBB = kruger_vars(X, A, B)
            
            dAdt = A*(rA*(A + RA*PBA) + B - RAA*PAA) - RA*PBA*(RAA*PAA + RB*PAB)
            dBdt = B*(rB*(B+RB*PAB) + A - RBB*PBB) - RB*PAB*(RBB*PBB + RA*PBA)
            FA = dAdt / (dAdt + dBdt)
            return [(fA - FA)/(1-X)]
        
        sol = solve_ivp(kruger_eqns, self.bounds, [self.fA0], method=self.method, max_step=0.001, t_eval=self.t_eval)
        
        X = sol.t
        fA = sol.y[0]
        
        A = fA*(1-X)*(self.fA0 + self.fB0)
        B = (1 - fA)*(1-X)*(self.fA0 + self.fB0)
        
        xA = (self.fA0 - A)/(self.fA0)
        xB = (self.fB0 - B)/(self.fB0)
        
        return X, xA, xB
    
    
def cpe_model(X_eval, r1, r2, b1, g1, b2, g2):
    
    f1 = 0.25
    inputs = [f1, r1, r2, b1, g1, b2, g2]
    copolymer = CopolymerEquations.fromLundberg(inputs)
    X1, xA1, xB1 = copolymer.solve_izu(X_eval)
    
    f1 = 0.5
    inputs = [f1, r1, r2, b1, g1, b2, g2]
    copolymer = CopolymerEquations.fromLundberg(inputs)
    X2, xA2, xB2 = copolymer.solve_izu(X_eval)
    
    f1 = 0.75
    inputs = [f1, r1, r2, b1, g1, b2, g2]
    copolymer = CopolymerEquations.fromLundberg(inputs)
    X3, xA3, xB3 = copolymer.solve_izu(X_eval)

    fit_data = np.concatenate((xA1, xB1, xA2, xB2, xA3, xB3))
    return fit_data

def cpe_model_(rng, r1, r2, b1, g1, b2, g2, size=None):
    
    X_eval = np.linspace(0.1, 0.7, 7)
    f1 = 0.25
    inputs = [f1, r1, r2, b1, g1, b2, g2]
    copolymer = CopolymerEquations.fromLundberg(inputs)

    # X_eval = np.linspace(0.1, 0.7, 7)
    X1, xA1, xB1 = copolymer.solve_izu(X_eval)

    fit_data = np.concatenate((xA1, xB1))
    return fit_data

# Modified function factory

from typing import Dict, Tuple, Optional

def create_wrapper(model_func, 
    params_fixed: Dict[str, float], 
    params_fit_guess: Dict[str, float], 
    ):
    
    fit_params = [p for p in params_fit_guess if p not in params_fixed]
    fit_guesses = [params_fit_guess[p] for p in fit_params]
    
    def fit_func(x, *args):
        params = {key: value for key, value in zip(fit_params, args)}
        params.update(params_fixed)
        #print(params)
        return model_func(x, **params)
    
    return fit_func

def cpe_irreversible_model_nll(log_r):
    
    log_r1, log_r2 = log_r
    r1 = 10.**float(log_r1)
    r2 = 10.**float(log_r2)
    
    
    r1_true = 3.0
    r2_true = 0.5
    
    b1_true = 0.0
    g1_true = 0.0
    b2_true = 0.0
    g2_true = 0.0
    
    X_eval = np.linspace(0.1, 0.7, 7)
    
    exp_data = cpe_model(X_eval, r1_true, r2_true, b1_true, g1_true, b2_true, g2_true)
    fit_data = cpe_model(X_eval, r1, r2, b1_true, g1_true, b2_true, g2_true)
    
    ssr = np.sum((exp_data - fit_data)**2)
    return np.log(ssr)
    # nll_ssr = -np.log(ssr)
    # return nll_ssr
    
    
from PyPESTO.FRP import create_FRP2_v1

def amici_CPE_sim(model, X_eval, f1, r1, r2, b1, g1, b2, g2):
    
    model.setParameterByName('kpAA', 1)
    model.setParameterByName('rA', r1)
    model.setParameterByName('rB', r2)
    model.setParameterByName('rX', 1)

    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-12)
    model.setInitialStates([0.000, 0.005, f1, 1-f1, 0, 0,0, 0, 0,  0])
    model.setTimepoints(np.linspace(0, 1000, 1000))
    rdata = amici.runAmiciSimulation(model, solver)
    
    meas_a = rdata.by_id('A')
    meas_b = rdata.by_id('B')
    
    conv_a = (f1 - meas_a) / f1
    conv_b = (1 - f1 - meas_b) / (1 - f1)
    conv = 1 - meas_a - meas_b
    
    xA = np.interp(X_eval, conv, conv_a)
    xB = np.interp(X_eval, conv, conv_b)
    
    return X_eval, xA, xB

def amici_CPErev_sim(model, X_eval, f1, r1, r2, b1, g1, b2, g2):
    
    model.setParameterByName('kpAA', 1)
    model.setParameterByName('rA', r1)
    model.setParameterByName('rB', r2)
    model.setParameterByName('rX', 1.)
    

def amici_CPE_sim_(model, X_eval, f1, r1, r2, b1, g1, b2, g2, rX=1.):
    
    model.setParameterByName('kpAA', 1)
    model.setParameterByName('rA', r1)
    model.setParameterByName('rB', r2)
    model.setParameterByName('rX', rX)

    solver = model.getSolver()
    solver.setAbsoluteTolerance(1e-12)
    model.setInitialStates([0.000, 0.005, f1, 1-f1, 0, 0,0, 0, 0,  0])
    model.setTimepoints(np.linspace(0, 5000, 5000))
    rdata = amici.runAmiciSimulation(model, solver)
    
    return rdata
    
def amici_CPE_model(X_eval, r1, r2, b1, g1, b2, g2):
    
    model_module = amici.import_model_module('FRP2_v1', 'tmp/FRP2_v1')
    model = model_module.getModel()

    f1 = 0.25
    X1, xA1, xB1 = amici_CPE_sim(model, X_eval, f1, r1, r2, b1, g1, b2, g2)
    fit_data = np.concatenate((xA1, xB1))
    return fit_data
    
    f1 = 0.5
    X2, xA2, xB2 = amici_CPE_sim(model, X_eval, f1, r1, r2, b1, g1, b2, g2)
    
    f1 = 0.75
    X3, xA3, xB3 = amici_CPE_sim(model, X_eval, f1, r1, r2, b1, g1, b2, g2)

    fit_data = np.concatenate((xA1, xB1, xA2, xB2, xA3, xB3))
    return fit_data

def amici_CPE_model(X_eval, r1, r2, b1, g1, b2, g2):
    
    model_module = amici.import_model_module('FRP2_v1', 'tmp/FRP2_v1')
    model = model_module.getModel()

    f1 = 0.25
    X1, xA1, xB1 = amici_CPE_sim(model, X_eval, f1, r1, r2, b1, g1, b2, g2)
    fit_data = np.concatenate((xA1, xB1))
    return fit_data


from functools import cache

# @cache
def amici_irr_acqf():
    
    r1_true = 4.0
    r2_true = 0.5
    
    b1_true = 0.0
    g1_true = 0.0
    b2_true = 0.0
    g2_true = 0.0
    
    X_eval = np.linspace(0.1, 0.7, 7)
    
    exp_data = amici_CPE_model(X_eval, r1_true, r2_true, b1_true, g1_true, b2_true, g2_true)
    return exp_data

def amici_irreversible_model(log_r):
    
    log_r1, log_r2 = log_r
    r1 = 10.**float(log_r1)
    r2 = 10.**float(log_r2)
    
    b1_true = 0.0
    g1_true = 0.0
    b2_true = 0.0
    g2_true = 0.0
    
    X_eval = np.linspace(0.1, 0.7, 7)
    
    exp_data = amici_irr_acqf()
    fit_data = amici_CPE_model(X_eval, r1, r2, b1_true, g1_true, b2_true, g2_true)
    
    ssr = np.sum((exp_data - fit_data)**2)
    return np.log(ssr)
    
    # create_FRP2_v1.run_amici_simulation(model, timepoints, cA0=f)