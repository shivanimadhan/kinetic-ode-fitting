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

        #############################################
        ###### Chain model (1st order moments) ######
        #############################################

        PAA_1 = y[11]
        PAB_1 = y[12]
        PBA_1 = y[13]
        PBB_1 = y[14]
        PA_1  = PAA_1 + PBA_1
        PB_1  = PAB_1 + PBB_1

        # 11 - PAA (1st order moment)
        dydt[11]  = self.kpAA*RA*A - self.kpAA*PAA_1*A - self.kpAB*PAA_1*B + \
                   self.kpAA*(PAA_0+PAA_1)*A + self.kpAA*(PBA_0+PBA_1)*A - \
                   self.kdAA*PAA_1 + self.kdAA*pPAAA*(PAA_1-PAA_0) + self.kdAB*pPAAB*(PAB_1-PAB_0) - \
                   self.ktAA*PAA_1*PA_0 - self.ktAB*PAA_1*PB_0

        # 12 - PAB (1st order moment)
        dydt[12] = self.kpAB*RA*B - self.kpBA*PAB_1*A - self.kpBB*PAB_1*B + \
                   self.kpAB*(PAA_0+PAA_1)*B + self.kpAB*(PBA_0+PBA_1)*B - \
                   self.kdAB*PAB_1 + self.kdBA*pPABA*(PBA_1-PBA_0) + self.kdBB*pPABB*(PBB_1-PBB_0) - \
                   self.ktAB*PAB_1*PA_0 - self.ktBB*PAB_1*PB_0

        # 13 - PBA (1st order moment)
        dydt[13] = self.kpBA*RB*A - self.kpAA*PBA_1*A - self.kpAB*PBA_1*B + \
                   self.kpBA*(PAB_0+PAB_1)*A + self.kpBA*(PBB_0+PBB_1)*A - \
                   self.kdBA*PBA_1 + self.kdAA*pPBAA*(PAA_1-PAA_0) + self.kdAB*pPBAB*(PAB_1-PAB_0) - \
                   self.ktAA*PBA_1*PA_0 - self.ktAB*PBA_1*PB_0

        # 14 - PBB (1st order moment)
        dydt[14] = self.kpBB*RB*B - self.kpBA*PBB_1*A - self.kpBB*PBB_1*B + \
                   self.kpBB*(PAB_0+PAB_1)*B + self.kpBB*(PBB_0+PBB_1)*B - \
                   self.kdBB*PBB_1 + self.kdBA*pPBBA*(PBA_1-PBA_0) + self.kdBB*pPBBB*(PBB_1-PBB_0) - \
                   self.ktAB*PBB_1*PA_0 - self.ktBB*PBB_1*PB_0
        
        # 15 - PD (1st order moment)
        dydt[15] = 0.5*self.ktcAA*(PA_0*PA_1+PA_1*PA_0) + 0.5*self.ktcAB*(PA_0*PB_1+PA_1*PB_0) + 0.5*self.ktcBB*(PB_0*PB_1+PB_1*PB_0) + \
                   self.ktdAA*PA_0*PA_1 + 0.5*self.ktdAB*(PA_0*PB_1+PB_0*PA_1) + self.ktdBB*PB_0*PB_1

        
        #############################################
        ###### Chain model (2nd order moments) ######
        #############################################

        PAA_2 = y[16]
        PAB_2 = y[17]
        PBA_2 = y[18]
        PBB_2 = y[19]
        PA_2 = PAA_2 + PBA_2
        PB_2 = PAB_2 + PBB_2

        # 16 - PAA (2nd order moment)
        dydt[16] = self.kpAA*RA*A - self.kpAA*PAA_2*A - self.kpAB*PAA_2*B + \
                   self.kpAA*(PAA_0+2*PAA_1+PAA_2)*A + self.kpAA*(PBA_0+2*PBA_1+PBA_2)*A - \
                   self.kdAA*PAA_2 + self.kdAA*pPAAA*(PAA_2-2*PAA_1+PAA_0) + self.kdAB*pPAAB*(PAB_2-2*PAB_1+PAB_0) - \
                   self.ktAA*PAA_2*PA_0 - self.ktAB*PAA_2*PB_0

        # 17 - PAB (2nd order moment)
        dydt[17] = self.kpAB*RA*B - self.kpBA*PAB_2*A - self.kpBB*PAB_2*B + \
                   self.kpAB*(PAA_0+2*PAA_1+PAA_2)*B + self.kpAB*(PBA_0+2*PBA_1+PBA_2)*B - \
                   self.kdAB*PAB_2 + self.kdBA*pPABA*(PBA_2-2*PBA_1+PBA_0) + self.kdBB*pPABB*(PBB_2-2*PBB_1+PBB_0) - \
                   self.ktAB*PAB_2*PA_0 - self.ktBB*PAB_2*PB_0

        # 18 - PBA (2nd order moment)
        dydt[18] = self.kpBA*RB*A - self.kpAA*PBA_2*A - self.kpAB*PBA_2*B + \
                   self.kpBA*(PAB_0+2*PAB_1+PAB_2)*A + self.kpBA*(PBB_0+2*PBB_1+PBB_2)*A - \
                   self.kdBA*PBA_2 + self.kdAA*pPBAA*(PAA_2-2*PAA_1+PAA_0) + self.kdAB*pPBAB*(PAB_2-2*PAB_1+PAB_0) - \
                   self.ktAA*PBA_2*PA_0 - self.ktAB*PBA_2*PB_0

        # 19 - PBB (2nd order moment)
        dydt[19] = self.kpBB*RB*B - self.kpBA*PBB_2*A - self.kpBB*PBB_2*B + \
                   self.kpBB*(PAB_0+2*PAB_1+PAB_2)*B + self.kpBB*(PBB_0+2*PBB_1+PBB_2)*B - \
                   self.kdBB*PBB_2 + self.kdBA*pPBBA*(PBA_2-2*PBA_1+PBA_0) + self.kdBB*pPBBB*(PBB_2-2*PBB_1+PBB_0) - \
                   self.ktAB*PBB_2*PA_0 - self.ktBB*PBB_2*PB_0
        
        # 20 - PD (2nd order moment)
        dydt[20] = 0.5*self.ktcAA*(PA_0*PA_2+2*PA_1*PA_1+PA_2*PA_0) + 0.5*self.ktcAB*(PA_0*PB_2+2*PA_1*PB_1+PA_2*PB_0) + 0.5*self.ktcBB*(PB_0*PB_2+2*PB_1*PB_1+PB_2*PB_0) + \
                   self.ktdAA*PA_0*PA_2 + 0.5*self.ktdAB*(PA_0*PB_2+PB_0*PA_2) + self.ktdBB*PB_0*PB_2

        ########################################################
        ###### Sequence model, active (0th order moments) ######
        ########################################################

        SA_0 = y[21]
        SB_0 = y[22]

        # 21 - Active A sequence (0th order moment)
        dydt[21] = self.kpAA*A*(R + SA_0 - SA_0) - self.kpAB*B*(RA + SA_0) + self.kpBA*(RB + SB_0)*A + \
                   self.kdAA*fPAA*SA_0 + self.kdAB*fPAB*SB_0 - self.kdAA*fPAA*SA_0 - self.kdBA*fPBA*SA_0 - \
                   self.ktAA*SA_0*SA_0 - self.ktAB*SA_0*SB_0

        # 22 - Active B sequence (0th order moment)
        dydt[22] = self.kpBB*B*(R + SB_0 - SB_0) - self.kpBA*A*(RB + SB_0) + self.kpAB*(RA + SA_0)*B + \
                   self.kdBB*fPBB*SB_0 + self.kdBA*fPBA*SA_0 - self.kdBB*fPBB*SB_0 - self.kdAB*fPAB*SB_0 - \
                   self.ktAB*SA_0*SB_0 - self.ktBB*SB_0*SB_0

        ########################################################
        ###### Sequence model, active (1st order moments) ######
        ########################################################

        SA_1 = y[23]
        SB_1 = y[24]

        # 23 - Active A sequence (1st order moment)
        dydt[23] = self.kpAA*A*(R + SA_0 + SA_1 - SA_1) - self.kpAB*B*(RA + SA_1) + self.kpBA*(RB + SB_0)*A + \
                   self.kdAA*fPAA*(SA_1-SA_0) + self.kdAB*fPAB*SB_0 - self.kdAA*fPAA*SA_1 - self.kdBA*fPBA*SA_0 - \
                   self.ktAA*SA_1*SA_0 - self.ktAB*SA_1*SB_0

        # 24 - Active B sequence (1st order moment)
        dydt[24] = self.kpBB*B*(R + SB_0 + SB_1 - SB_1) - self.kpBA*A*(RB + SB_1) + self.kpAB*(RA + SA_0)*B + \
                   self.kdBB*fPBB*(SB_1-SB_0) + self.kdBA*fPBA*SA_0 - self.kdBB*fPBB*SB_1 - self.kdAB*fPAB*SB_0 - \
                   self.ktAB*SA_0*SB_1 - self.ktBB*SB_0*SB_1
        
        ########################################################
        ###### Sequence model, active (2nd order moments) ######
        ########################################################

        SA_2 = y[25]
        SB_2 = y[26]

        # 25 - Active A sequence (2nd order moment)
        dydt[25] = self.kpAA*A*(R + SA_0 + 2*SA_1 + SA_2 - SA_2) - self.kpAB*B*(RA + SA_2) + self.kpBA*(RB + SB_0)*A + \
                   self.kdAA*fPAA*(SA_2 - 2*SA_1 + SA_0) + self.kdAB*fPAB*SB_0 - self.kdAA*fPAA*SA_2 - self.kdBA*fPBA*SA_0 - \
                   self.ktAA*SA_2*SA_0 - self.ktAB*SA_2*SB_0

        # 26 - Active B sequence (2nd order moment)
        dydt[26] = self.kpBB*B*(R + SB_0 + 2*SB_1 + SB_2 - SB_2) - self.kpBA*A*(RB + SB_2) + self.kpAB*(RA + SA_0)*B + \
                   self.kdBB*fPBB*(SB_2 - 2*SB_1 + SB_0) + self.kdBA*fPBA*SA_0 - self.kdBB*fPBB*SB_2 - self.kdAB*fPAB*SB_0 - \
                   self.ktAB*SA_0*SB_2 - self.ktBB*SB_0*SB_2
        

        ##########################################################
        ###### Sequence model, inactive (0th order moments) ######
        ##########################################################

        # 27 - Inactive A sequence (0th order moment)
        dydt[27] = self.kpAB*B*(RA + SA_0) - self.kdAB*fPAB*SB_0 + \
                   self.ktAA*SA_0*SA_0 + self.ktAB*SA_0*SB_0

        # 28 - Inactive B sequence (0th order moment)
        dydt[28] = self.kpBA*A*(RB + SB_0) - self.kdBA*fPBA*SA_0 + \
                   self.ktAB*SA_0*SB_0 + self.ktBB*SB_0*SB_0

        ##########################################################
        ###### Sequence model, inactive (1st order moments) ######
        ##########################################################

        # 29 - Inactive A sequence (1st order moment)
        dydt[29] = self.kpAB*B*(RA + SA_1) - self.kdAB*fPAB*SB_0 + \
                   self.ktcAA*(2*SA_0*SA_1) + self.ktdAA*SA_0*SA_1 + self.ktAB*SA_1*SB_0

        # 30 - Inactive B sequence (1st order moment)
        dydt[30] = self.kpBA*A*(RB + SB_1) - self.kdBA*fPBA*SA_0 + \
                   self.ktAB*SA_0*SB_1 + self.ktcBB*(2*SB_0*SB_1) + self.ktdBB*SB_0*SB_1

        ##########################################################
        ###### Sequence model, inactive (2nd order moments) ######
        ##########################################################

        # 31 - Inactive A sequence (2nd order moment)
        dydt[31] = self.kpAB*B*(RA + SA_2) - self.kdAB*fPAB*SB_0 + \
                   self.ktcAA*(2*SA_0*SA_2+2*SA_1*SA_1) + self.ktdAA*SA_0*SA_2 + self.ktAB*SA_2*SB_0


        # 32 - Inactive B sequence (2nd order moment)
        dydt[32] = self.kpBA*A*(RB + SB_2) - self.kdBA*fPBA*SA_0 + \
                   self.ktAB*SA_0*SB_2 + self.ktcBB*(2*SB_0*SB_2+2*SB_0*SB_1) + self.ktdBB*SB_0*SB_2

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

        # Define monomer concentrations
        self.A = self.sol[2] - 1e-10
        self.B = self.sol[3] - 1e-10

        # Define copolymer composition
        self.CA = (self.A0 - self.A) / (self.A0 - self.A + self.B0 - self.B)
        self.CB = (self.B0 - self.B) / (self.A0 - self.A + self.B0 - self.B)

        # Define chain model moments
        self.P_0 = self.sol[6]  + self.sol[7]  + self.sol[8]  + self.sol[9]  + self.sol[10] + 1e-40
        self.P_1 = self.sol[11] + self.sol[12] + self.sol[13] + self.sol[14] + self.sol[15] + 1e-40
        self.P_2 = self.sol[16] + self.sol[17] + self.sol[18] + self.sol[19] + self.sol[20] + 1e-40
        self.PD_0 = self.sol[10] + 1e-40
        self.PD_1 = self.sol[15] + 1e-40
        self.PD_2 = self.sol[20] + 1e-40


        # Calculate conversion
        self.x_0  = self.P_1 / (self.P_1 + self.A + self.B) # Chain model
        self.x_1  = (self.A0 + self.B0 - self.A - self.B) / (self.A0 + self.B0) # Concentration
        self.x    = self.x_1
        self.xA   = (self.A0 - self.A) / self.A0
        self.xB   = (self.B0 - self.B) / self.B0

        # Calculate average chain model parameters
        self.NACL = self.P_1 / self.P_0
        self.WACL = self.P_2 / self.P_1
        self.PDI  = self.WACL / self.NACL

        # Define sequence model moments
        self.SA_0 = self.sol[21] + self.sol[27] + 1e-40
        self.SB_0 = self.sol[22] + self.sol[28] + 1e-40
        self.SA_1 = self.sol[23] + self.sol[29] + 1e-40
        self.SB_1 = self.sol[24] + self.sol[30] + 1e-40
        self.SA_2 = self.sol[25] + self.sol[31] + 1e-40
        self.SB_2 = self.sol[26] + self.sol[32] + 1e-40

        # Calculate average sequence model parameters
        self.NASL_A = self.SA_1 / self.SA_0
        self.NASL_B = self.SB_1 / self.SB_0
        self.WASL_A = self.SA_2 / self.SA_1
        self.WASL_B = self.SB_2 / self.SB_1
        self.SDI_A  = self.WASL_A / self.NASL_A
        self.SDI_B  = self.WASL_B / self.NASL_B
        
        
