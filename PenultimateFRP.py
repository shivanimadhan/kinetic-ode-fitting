import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

class CopolymerizationModel():

    def method_of_moments_ODEs(self, t, y):
        dydt = np.zeros_like(y)

        # Define chain end dyad concentrations
        PAA_0 = y[6]
        PAB_0 = y[7]
        PBA_0 = y[8]
        PBB_0 = y[9]
        PA_0 = PAA_0 + PBA_0 + 1e-40
        PB_0 = PAB_0 + PBB_0 + 1e-40

        # Calculate chian end dyad fractions
        fPAA = PAA_0 / PA_0
        fPAB = PAB_0 / PB_0
        fPBA = PBA_0 / PA_0
        fPBB = PBB_0 / PB_0
        
        # Define chain end triad fractions from Markov chain factorization
        pPAAA = pPAAB = fPAA
        pPBAA = pPBAB = fPBA
        pPABA = pPABB = fPAB
        pPBBA = pPBBB = fPBB

        ###########################################
        ###### Small molecule concentrations ######
        ###########################################

        I  = y[0]
        R  = y[1]
        A  = y[2]
        B  = y[3]
        RA = y[4]
        RB = y[5]

        x  = (self.A0 + self.B0 - A - B) / (self.A0 + self.B0)

        if x > (1 - 1e-3):
            return dydt

        # 00 - I, Initiator concentration
        dydt[0]  = -self.kd*I

        # 01 - R, Radical concentration
        dydt[1]  = 2*self.f*self.kd*I - self.kpRA*R*A - self.kpRB*R*B

        # 02 - A monomer concentration
        dydt[2]  = -A*(self.kpRA*R + self.kpAA*RA + self.kpBA*RB + self.kpAA*PA_0 + self.kpBA*PB_0) + self.kdAA*PAA_0 + self.kdBA*PBA_0

        # 03 - B monomer concentration
        dydt[3]  = -B*(self.kpRB*R + self.kpAB*RA + self.kpBB*RB + self.kpAB*PA_0 + self.kpBB*PB_0) + self.kdAB*PAB_0 + self.kdBB*PBB_0

        # 04 - RA dyad concentration
        dydt[4]  = self.kpRA*R*A - self.kpAA*RA*A - self.kpAB*RA*B

        # 05 - RB dyad concentration
        dydt[5]  = self.kpRB*R*B - self.kpBB*RB*B - self.kpBA*RB*A

        #############################################
        ###### Chain model (0th order moments) ######
        #############################################

        PD_0 = y[10]

        # 06 - PAA (0th order moment) 
        dydt[6]  = self.kpAA*RA*A - self.kpAA*PAA_0*A - self.kpAB*PAA_0*B + \
                   self.kpAA*PAA_0*A + self.kpAA*PBA_0*A - \
                   self.kdAA*PAA_0 + self.kdAA*pPAAA*PAA_0 + self.kdAB*pPAAB*PAB_0 - \
                   self.ktAA*PAA_0*PA_0 - self.ktAB*PAA_0*PB_0

        # 07 - PAB (0th order moment)
        dydt[7]  = self.kpAB*RA*B - self.kpBA*PAB_0*A - self.kpBB*PAB_0*B + \
                   self.kpAB*PAA_0*B + self.kpAB*PBA_0*B - \
                   self.kdAB*PAB_0 + self.kdBA*pPABA*PBA_0 + self.kdBB*pPABB*PBB_0 - \
                   self.ktAB*PAB_0*PAA_0 - self.ktBB*PAB_0*PB_0

        # 08 - PBA (0th order moment)
        dydt[8]  = self.kpBA*RB*A - self.kpAA*PBA_0*A - self.kpAB*PBA_0*B + \
                   self.kpBA*PAB_0*A + self.kpBA*PBB_0*A - \
                   self.kdBA*PBA_0 + self.kdAA*pPBAA*PAA_0 + self.kdAB*pPBAB*PAB_0 - \
                   self.ktAA*PBA_0*PA_0 - self.ktAB*PBA_0*PB_0

        # 09 - PBB (0th order moment)
        dydt[9]  = self.kpBB*RB*B - self.kpBA*PBB_0*A - self.kpBB*PBB_0*B + \
                   self.kpBB*PAB_0*B + self.kpBB*PBB_0*B - \
                   self.kdBB*PBB_0 + self.kdBA*pPBBA*PBA_0 + self.kdBB*pPBBB*PBB_0 - \
                   self.ktAB*PBB_0*PA_0 - self.ktBB*PBB_0*PB_0
        
        # 10 - PD (0th order moment)
        dydt[10] = 0.5*self.ktcAA*PA_0*PA_0 + 0.5*self.ktcAB*PA_0*PB_0 + 0.5*self.ktcBB*PB_0*PB_0 + \
                   self.ktdAA*PA_0*PA_0 + 0.5*self.ktdAB*(PA_0*PB_0 + PB_0*PA_0) + self.ktdBB*PB_0*PB_0

        return dydt

    def __init__(self, k, y0, t_span):

        # Define kinetic rate constants
        (
            self.kd, self.f, # Initiator dissociation rate constant
            self.kpAA, self.kpAB, self.kpBA, self.kpBB, # Propagation rate constants
            self.kdAA, self.kdAB, self.kdBA, self.kdBB, # Depropagation rate constants
            self.ktcAA, self.ktcAB, self.ktcBB, # Termination (combination) rate constants
            self.ktdAA, self.ktdAB, self.ktdBB  # Termination (disproportionation) rate constants
        ) = k
        self.kpRA = self.kpAA
        self.kpRB = self.kpBB
        self.ktAA = self.ktcAA + self.ktdAA
        self.ktAB = self.ktcAB + self.ktdAB
        self.ktBB = self.ktcBB + self.ktdBB

        # Define initial values
        self.y0 = y0
        self.I0 = self.y0[0]
        self.A0 = self.y0[2]
        self.B0 = self.y0[3]
        
        # Solve the method of moments equations
        self.t_span = t_span
        method = None
        if self.kdAA*self.kdAB*self.kdBA*self.kdBB < 1e-5:
            method = 'BDF'
        else:
            method = 'LSODA'
        self.mom_sol = solve_ivp(self.method_of_moments_ODEs, self.t_span, self.y0, method=method, rtol=1e-10, atol=1e-10, first_step=1e-10)

        # Interpolate method of moments solution
        num_points = 40
        #self.t = np.linspace(np.min(self.mom_sol.t), np.max(self.mom_sol.t), num_points)t_span = [0, 60*3600]
        self.t = np.linspace(self.t_span[0], self.t_span[1], num_points)
        self.t_hours = self.t / 3600
        self.sol = np.zeros([len(self.mom_sol.y), num_points])

        for i in range(len(self.sol)):
            interp_func = interp1d(self.mom_sol.t, self.mom_sol.y[i], bounds_error=False, kind='cubic')
            self.sol[i,:] = interp_func(self.t)

        self.sol[np.isnan(self.sol)] = 0

        # Define monomer concentrations
        self.A = self.sol[2]
        self.B = self.sol[3] 

        # Define copolymer composition
        self.CA = (self.A0 - self.A) / (self.A0 - self.A + self.B0 - self.B - 1e-40)
        self.CB = (self.B0 - self.B) / (self.A0 - self.A + self.B0 - self.B - 1e-40)

        # Define chain model moments
        self.P_0 = self.sol[6]  + self.sol[7]  + self.sol[8]  + self.sol[9]  + self.sol[10] + 1e-40
        self.PD_0 = self.sol[10] + 1e-40


        # Calculate conversion
        self.x_1  = (self.A0 + self.B0 - self.A - self.B) / (self.A0 + self.B0) # Concentration
        self.x    = self.x_1
        self.xA   = (self.A0 - self.A) / self.A0
        self.xB   = (self.B0 - self.B) / self.B0
        
        
