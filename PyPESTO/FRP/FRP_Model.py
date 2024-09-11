import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from typing import List, Optional, Tuple

class CopolymerizationModel():

    def kinetic_odes(self, t, y):
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
        
        x = (self.A0 + self.B0 - A - B) / (self.A0 + self.B0)
        # print(x)
        if (x > (1 - 1e-4)):
            # print(x)
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
                   2*self.ktAA*PAA_0*(PAA_0+1.0*PBA_0) - self.ktAB*PAA_0*PB_0

        # 07 - PAB (0th order moment)
        dydt[7]  = self.kpAB*RA*B - self.kpBA*PAB_0*A - self.kpBB*PAB_0*B + \
                   self.kpAB*PAA_0*B + self.kpAB*PBA_0*B - \
                   self.kdAB*PAB_0 + self.kdBA*pPABA*PBA_0 + self.kdBB*pPABB*PBB_0 - \
                   self.ktAB*PAB_0*PA_0 - 2*self.ktBB*PAB_0*(PAB_0+1.0*PBB_0)

        # 08 - PBA (0th order moment)
        dydt[8]  = self.kpBA*RB*A - self.kpAA*PBA_0*A - self.kpAB*PBA_0*B + \
                   self.kpBA*PAB_0*A + self.kpBA*PBB_0*A - \
                   self.kdBA*PBA_0 + self.kdAA*pPBAA*PAA_0 + self.kdAB*pPBAB*PAB_0 - \
                   2*self.ktAA*PBA_0*(1.0*PAA_0+PBA_0) - self.ktAB*PBA_0*PB_0

        # 09 - PBB (0th order moment)
        dydt[9]  = self.kpBB*RB*B - self.kpBA*PBB_0*A - self.kpBB*PBB_0*B + \
                   self.kpBB*PAB_0*B + self.kpBB*PBB_0*B - \
                   self.kdBB*PBB_0 + self.kdBA*pPBBA*PBA_0 + self.kdBB*pPBBB*PBB_0 - \
                   self.ktAB*PBB_0*PA_0 - 2*self.ktBB*PBB_0*(1.0*PAB_0+PBB_0)
        
        # 10 - PD (0th order moment)
        dydt[10] = self.ktcAA*PA_0*PA_0 + self.ktcAB*PA_0*PB_0 + self.ktcBB*PB_0*PB_0 + \
                   2*self.ktdAA*PA_0*PA_0 + 2*0.5*self.ktdAB*(PA_0*PB_0 + PB_0*PA_0) + 2*self.ktdBB*PB_0*PB_0

        #############################################
        ###### Chain model (1st order moments) ######
        #############################################

        if not self.solve_chain_model:
            return dydt

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
                   2*self.ktAA*PAA_1*PA_0 - self.ktAB*PAA_1*PB_0

        # 12 - PAB (1st order moment)
        dydt[12] = self.kpAB*RA*B - self.kpBA*PAB_1*A - self.kpBB*PAB_1*B + \
                   self.kpAB*(PAA_0+PAA_1)*B + self.kpAB*(PBA_0+PBA_1)*B - \
                   self.kdAB*PAB_1 + self.kdBA*pPABA*(PBA_1-PBA_0) + self.kdBB*pPABB*(PBB_1-PBB_0) - \
                   self.ktAB*PAB_1*PA_0 - 2*self.ktBB*PAB_1*PB_0

        # 13 - PBA (1st order moment)
        dydt[13] = self.kpBA*RB*A - self.kpAA*PBA_1*A - self.kpAB*PBA_1*B + \
                   self.kpBA*(PAB_0+PAB_1)*A + self.kpBA*(PBB_0+PBB_1)*A - \
                   self.kdBA*PBA_1 + self.kdAA*pPBAA*(PAA_1-PAA_0) + self.kdAB*pPBAB*(PAB_1-PAB_0) - \
                   2*self.ktAA*PBA_1*PA_0 - self.ktAB*PBA_1*PB_0

        # 14 - PBB (1st order moment)
        dydt[14] = self.kpBB*RB*B - self.kpBA*PBB_1*A - self.kpBB*PBB_1*B + \
                   self.kpBB*(PAB_0+PAB_1)*B + self.kpBB*(PBB_0+PBB_1)*B - \
                   self.kdBB*PBB_1 + self.kdBA*pPBBA*(PBA_1-PBA_0) + self.kdBB*pPBBB*(PBB_1-PBB_0) - \
                   self.ktAB*PBB_1*PA_0 - 2*self.ktBB*PBB_1*PB_0
        
        # 15 - PD (1st order moment)
        dydt[15] = self.ktcAA*(PA_0*PA_1+PA_1*PA_0) + self.ktcAB*(PA_0*PB_1+PA_1*PB_0) + self.ktcBB*(PB_0*PB_1+PB_1*PB_0) + \
                   2*self.ktdAA*PA_0*PA_1 + 2*0.5*self.ktdAB*(PA_0*PB_1+PB_0*PA_1) + 2*self.ktdBB*PB_0*PB_1

        
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
                   2*self.ktAA*PAA_2*PA_0 - self.ktAB*PAA_2*PB_0

        # 17 - PAB (2nd order moment)
        dydt[17] = self.kpAB*RA*B - self.kpBA*PAB_2*A - self.kpBB*PAB_2*B + \
                   self.kpAB*(PAA_0+2*PAA_1+PAA_2)*B + self.kpAB*(PBA_0+2*PBA_1+PBA_2)*B - \
                   self.kdAB*PAB_2 + self.kdBA*pPABA*(PBA_2-2*PBA_1+PBA_0) + self.kdBB*pPABB*(PBB_2-2*PBB_1+PBB_0) - \
                   self.ktAB*PAB_2*PA_0 - 2*self.ktBB*PAB_2*PB_0

        # 18 - PBA (2nd order moment)
        dydt[18] = self.kpBA*RB*A - self.kpAA*PBA_2*A - self.kpAB*PBA_2*B + \
                   self.kpBA*(PAB_0+2*PAB_1+PAB_2)*A + self.kpBA*(PBB_0+2*PBB_1+PBB_2)*A - \
                   self.kdBA*PBA_2 + self.kdAA*pPBAA*(PAA_2-2*PAA_1+PAA_0) + self.kdAB*pPBAB*(PAB_2-2*PAB_1+PAB_0) - \
                   2*self.ktAA*PBA_2*PA_0 - self.ktAB*PBA_2*PB_0

        # 19 - PBB (2nd order moment)
        dydt[19] = self.kpBB*RB*B - self.kpBA*PBB_2*A - self.kpBB*PBB_2*B + \
                   self.kpBB*(PAB_0+2*PAB_1+PAB_2)*B + self.kpBB*(PBB_0+2*PBB_1+PBB_2)*B - \
                   self.kdBB*PBB_2 + self.kdBA*pPBBA*(PBA_2-2*PBA_1+PBA_0) + self.kdBB*pPBBB*(PBB_2-2*PBB_1+PBB_0) - \
                   self.ktAB*PBB_2*PA_0 - 2*self.ktBB*PBB_2*PB_0
        
        # 20 - PD (2nd order moment)
        dydt[20] = self.ktcAA*(PA_0*PA_2+2*PA_1*PA_1+PA_2*PA_0) + self.ktcAB*(PA_0*PB_2+2*PA_1*PB_1+PA_2*PB_0) + self.ktcBB*(PB_0*PB_2+2*PB_1*PB_1+PB_2*PB_0) + \
                   2*self.ktdAA*PA_0*PA_2 + 2*0.5*self.ktdAB*(PA_0*PB_2+PB_0*PA_2) + 2*self.ktdBB*PB_0*PB_2

        ########################################################
        ###### Sequence model, active (0th order moments) ######
        ########################################################

        if not self.solve_sequence_model:
            return dydt

        SA_0 = y[21]
        SB_0 = y[22]

        # 21 - Active A sequence (0th order moment)
        dydt[21] = self.kpAA*A*(R + SA_0 - SA_0) - self.kpAB*B*(RA + SA_0) + self.kpBA*A*(RB + SB_0) + \
                   self.kdAA*fPAA*SA_0 + self.kdAB*fPAB*SB_0 - self.kdAA*fPAA*SA_0 - self.kdBA*fPBA*SA_0 - \
                   2*self.ktAA*SA_0*SA_0 - self.ktAB*SA_0*SB_0

        # 22 - Active B sequence (0th order moment)
        dydt[22] = self.kpBB*B*(R + SB_0 - SB_0) - self.kpBA*A*(RB + SB_0) + self.kpAB*(RA + SA_0)*B + \
                   self.kdBB*fPBB*SB_0 + self.kdBA*fPBA*SA_0 - self.kdBB*fPBB*SB_0 - self.kdAB*fPAB*SB_0 - \
                   self.ktAB*SA_0*SB_0 - 2*self.ktBB*SB_0*SB_0

        ########################################################
        ###### Sequence model, active (1st order moments) ######
        ########################################################

        SA_1 = y[23]
        SB_1 = y[24]

        # 23 - Active A sequence (1st order moment)
        dydt[23] = self.kpAA*A*(R + SA_0 + SA_1 - SA_1) - self.kpAB*B*(RA + SA_1) + self.kpBA*(RB + SB_0)*A + \
                   self.kdAA*fPAA*(SA_1-SA_0) + self.kdAB*fPAB*SB_0 - self.kdAA*fPAA*SA_1 - self.kdBA*fPBA*SA_0 - \
                   2*self.ktAA*SA_1*SA_0 - self.ktAB*SA_1*SB_0

        # 24 - Active B sequence (1st order moment)
        dydt[24] = self.kpBB*B*(R + SB_0 + SB_1 - SB_1) - self.kpBA*A*(RB + SB_1) + self.kpAB*(RA + SA_0)*B + \
                   self.kdBB*fPBB*(SB_1-SB_0) + self.kdBA*fPBA*SA_0 - self.kdBB*fPBB*SB_1 - self.kdAB*fPAB*SB_0 - \
                   self.ktAB*SA_0*SB_1 - 2*self.ktBB*SB_0*SB_1
        
        ########################################################
        ###### Sequence model, active (2nd order moments) ######
        ########################################################

        SA_2 = y[25]
        SB_2 = y[26]

        # 25 - Active A sequence (2nd order moment)
        dydt[25] = self.kpAA*A*(R + SA_0 + 2*SA_1 + SA_2 - SA_2) - self.kpAB*B*(RA + SA_2) + self.kpBA*(RB + SB_0)*A + \
                   self.kdAA*fPAA*(SA_2 - 2*SA_1 + 1*SA_0) + self.kdAB*fPAB*SB_0 - self.kdAA*fPAA*SA_2 - self.kdBA*fPBA*SA_0 - \
                   2*self.ktAA*SA_2*SA_0 - self.ktAB*SA_2*SB_0

        # 26 - Active B sequence (2nd order moment)
        dydt[26] = self.kpBB*B*(R + SB_0 + 2*SB_1 + SB_2 - SB_2) - self.kpBA*A*(RB + SB_2) + self.kpAB*(RA + SA_0)*B + \
                   self.kdBB*fPBB*(SB_2 - 2*SB_1 + 1*SB_0) + self.kdBA*fPBA*SA_0 - self.kdBB*fPBB*SB_2 - self.kdAB*fPAB*SB_0 - \
                   self.ktAB*SA_0*SB_2 - 2*self.ktBB*SB_0*SB_2
        

        ##########################################################
        ###### Sequence model, inactive (0th order moments) ######
        ##########################################################

        # 27 - Inactive A sequence (0th order moment)
        dydt[27] = self.kpAB*B*(RA + SA_0) - self.kdAB*fPAB*SB_0 + \
                   2*self.ktAA*SA_0*SA_0 + self.ktAB*SA_0*SB_0

        # 28 - Inactive B sequence (0th order moment)
        dydt[28] = self.kpBA*A*(RB + SB_0) - self.kdBA*fPBA*SA_0 + \
                   self.ktAB*SA_0*SB_0 + 2*self.ktBB*SB_0*SB_0

        ##########################################################
        ###### Sequence model, inactive (1st order moments) ######
        ##########################################################

        # 29 - Inactive A sequence (1st order moment)
        dydt[29] = self.kpAB*B*(RA + SA_1) - self.kdAB*fPAB*SB_0 + \
                   self.ktcAA*(2*SA_0*SA_1) + 2*self.ktdAA*SA_0*SA_1 + self.ktAB*SA_1*SB_0

        # 30 - Inactive B sequence (1st order moment)
        dydt[30] = self.kpBA*A*(RB + SB_1) - self.kdBA*fPBA*SA_0 + \
                   self.ktAB*SA_0*SB_1 + self.ktcBB*(2*SB_0*SB_1) + 2*self.ktdBB*SB_0*SB_1

        ##########################################################
        ###### Sequence model, inactive (2nd order moments) ######
        ##########################################################

        # 31 - Inactive A sequence (2nd order moment)
        dydt[31] = self.kpAB*B*(RA + SA_2) - self.kdAB*fPAB*SB_0 + \
                   self.ktcAA*(2*SA_0*SA_2+2*SA_1*SA_1) + 2*self.ktdAA*SA_0*SA_2 + self.ktAB*SA_2*SB_0


        # 32 - Inactive B sequence (2nd order moment)
        dydt[32] = self.kpBA*A*(RB + SB_2) - self.kdBA*fPBA*SA_0 + \
                   self.ktAB*SA_0*SB_2 + self.ktcBB*(2*SB_0*SB_2+2*SB_1*SB_1) + 2*self.ktdBB*SB_0*SB_2

        return dydt

    def get_rxn_strings(self): #-> Tuple[List[str], List[str]]:
        # Define all polymer species
        polymer_species = ["PR", "PRA", "PRB", "PAA", "PAB", "PBA", "PBB", "PD"]
        
        # 35 reactions
        reaction_strings = [
            # Initiator decomposition
            "AIBN -kd-> PR",

            # Propagation
            "PR + A -kpAA-> PRA",
            "PR + B -kpBB-> PRB",

            "PRA + A -kpAA-> PAA",
            "PRA + B -kpAB-> PAB",
            "PRB + A -kpBA-> PBA",
            "PRB + B -kpBB-> PBB",

            "PAA + A -kpAA-> PAA",
            "PAA + B -kpAB-> PAB",
            "PAB + A -kpBA-> PBA",
            "PAB + B -kpBB-> PBB",

            "PBA + A -kpAA-> PAA",
            "PBA + B -kpAB-> PAB",
            "PBB + A -kpBA-> PBA",
            "PBB + B -kpBB-> PBB",

            # Depropagation
            "PAA -kdAA-> PAA",
            "PAB -kdAB-> PAB",
            "PBA -kdBA-> PBA",
            "PBB -kdBB-> PBB",

            # Termination by Combination
            "PAA -ktcAA-> PD",
            "PBA -ktcAA-> PD",

            "PAA -ktcAB-> PD",
            "PBA -ktcAB-> PD",
            "PAB -ktcAB-> PD",
            "PBB -ktcAB-> PD",

            "PAB -ktcBB-> PD",
            "PBB -ktcBB-> PD",

            # Termination by Disproportionation
            "PAA -ktdAA-> PD",
            "PBA -ktdAA-> PD",

            "PAA -ktdAB-> PD",
            "PBA -ktdAB-> PD",
            "PAB -ktdAB-> PD",
            "PBB -ktdAB-> PD",

            "PAB -ktdBB-> PD",
            "PBB -ktdBB-> PD"
        ]
        
        # Remove initiator decomposition reaction
        reaction_strings = reaction_strings[1:len(reaction_strings)]

        return polymer_species, reaction_strings


    def save_rxn_rates(self, output_filepath:str='Probs.txt'):

        t = self.t
        y = self.sol

        r = np.zeros((36, len(y[0])))

        I   = y[0,:]
        R   = y[1,:]
        A   = y[2,:]
        B   = y[3,:]
        RA  = y[4,:]
        RB  = y[5,:]
        PAA = y[6,:]
        PAB = y[7,:]
        PBA = y[8,:]
        PBB = y[9,:]
        PA = PAA + PBA
        PB = PAB + PBB

        # f = 2.0
        # # self.ktcAA = f*self.ktcAA
        # self.ktdAA = f*self.ktdAA
        # # self.ktcAB = f/2*self.ktcAB
        # self.ktdAB = f*self.ktdAB
        # # self.ktcBB = f*self.ktcBB
        # self.ktdBB = f*self.ktdBB

        # Time points
        r[0,:] = t

        ######################
        ##### Initiation #####
        ######################
        
        # Initiator decomposition: I -> 2R
        r[1,:]  = self.kd*I

        # Initiation: R + A -> RA
        r[2,:]  = self.kpAA*R*A

        # Initiation: R + B -> RB
        r[3,:]  = self.kpBB*R*B


        #############################
        ##### First Propagation #####
        #############################

        # Propagation: RA + A -kpAA-> PAA
        r[4,:]  = self.kpAA*RA*A

        # Propagation: RA + B -kpAB-> PAB
        r[5,:]  = self.kpAB*RA*B

        # Propagation: RB + A -kpBA-> PBA
        r[6,:]  = self.kpBA*RB*A
         
        # Propagation: RB + B -kpBB-> PBB
        r[7,:]  = self.kpBB*RB*B


        #######################
        ##### Propagation #####
        #######################

        # Propagation: PAA + A -kpAA-> PAA
        r[8,:]  = self.kpAA*PAA*A

        # Propagation: PAA + B -kpAB-> PAB
        r[9,:]  = self.kpAB*PAA*B

        # Propagation: PAB + A -kpBA-> PBA
        r[10,:] = self.kpBA*PAB*A

        # Propagation: PAB + B -kpBB-> PBB
        r[11,:] = self.kpBB*PAB*B 

        # Propagation: PBA + A -kpAA-> PAA
        r[12,:] = self.kpAA*PBA*A

        # Propagation: PBA + B -kpAB-> PAB
        r[13,:] = self.kpAB*PBA*B 

        # Propagation: PBB + A -kpBA-> PBA
        r[14,:] = self.kpBA*PBB*A

        # Propagation: PBB + B -kpBB-> PBB
        r[15,:] = self.kpBB*PBB*B
        

        #########################
        ##### Depropagation #####
        #########################

        # Depropagation: PAA -kdAA-> PA + A
        r[16, :] = self.kdAA*PAA

        # Depropagation: PAB -kdAB-> PA + B
        r[17, :] = self.kdAB*PAB

        # Depropagation: PBA -kdBA-> PB + A
        r[18, :] = self.kdBA*PBA

        # Depropagation: PBB -kdBB-> PB + B
        r[19, :] = self.kdBB*PBB


        ######################################
        ##### Termination by Combination #####
        ######################################

        # Termination by comb.: PAA + PA -> D
        r[20,:] = self.ktcAA*PAA*PA

        # Termination by comb.: PBA + PA -> D
        r[21,:] = self.ktcAA*PBA*PA


        # Termination by comb.: PAA + PB -> D
        r[22,:] = self.ktcAB*PAA*PB

        # Termination by comb.: PBA + PB -> D
        r[23,:] = self.ktcAB*PBA*PB

        # Termination by comb.: PA + PAB -> D
        r[24,:] = self.ktcAB*PA*PAB

        # Termination by comb.: PA + PBB -> D
        r[25,:] = self.ktcAB*PA*PBB


        # Termination by comb.: PAB + PB -> D
        r[26,:] = self.ktcBB*PAB*PB

        # Termination by comb.: PBB + PB -> D
        r[27,:] = self.ktcBB*PBB*PB


        #############################################
        ##### Termination by Disproportionation #####
        #############################################
        
        f = 2.0

        # Termination by disp.: PAA + PA -> 2D
        r[28,:] = f*self.ktdAA*PAA*PA

        # Termination by disp.: PBA + PA -> 2D
        r[29,:] = f*self.ktdAA*PBA*PA


        # Termination by disp.: PAA + PB -> 2D
        r[30,:] = f*self.ktdAB*PAA*PB/2

        # Termination by disp.: PBA + PB -> 2D
        r[31,:] = f*self.ktdAB*PBA*PB/2

        # Termination by disp.: PA + PAB -> 2D
        r[32,:] = f*self.ktdAB*PA*PAB/2

        # Termination by disp.: PA + PBB -> 2D
        r[33,:] = f*self.ktdAB*PA*PBB/2


        # Termination by disp.: PAB + PB -> 2D
        r[34,:] = f*self.ktdBB*PAB*PB

        # Termination by disp.: PBB + PB -> 2D
        r[35,:] = f*self.ktdBB*PBB*PB
        
        # r[29:,:] = 0
        # r[3,:] = 0
        # r[5:8,:] = 0
        # r[9:16,:] = 0

        # # Termination by disp.: PAA + PA -> 2D
        # r[28,:] = np.sqrt(self.ktdAA*PAA*PA)

        # # Termination by disp.: PBA + PA -> 2D
        # r[29,:] = np.sqrt(self.ktdAA*PBA*PA)

        # # Termination by disp.: PAA + PB -> 2D
        # r[30,:] = np.sqrt(self.ktdAB*PAA*PB)

        # # Termination by disp.: PBA + PB -> 2D
        # r[31,:] = np.sqrt(self.ktdAB*PBA*PB)

        # # Termination by disp.: PA + PAB -> 2D
        # r[32,:] = np.sqrt(self.ktdAB*PA*PAB)

        # # Termination by disp.: PA + PBB -> 2D
        # r[33,:] = np.sqrt(self.ktdAB*PA*PBB)

        # # Termination by disp.: PAB + PB -> 2D
        # r[34,:] = np.sqrt(self.ktdBB*PAB*PB)

        # # Termination by disp.: PBB + PB -> 2D
        # r[35,:] = np.sqrt(self.ktdBB*PBB*PB)

        # self.ktcAA = 2*self.ktcAA
        # self.ktdAA = 2*self.ktdAA
        # self.ktcAB = 2*self.ktcAB
        # self.ktdAB = 2*self.ktdAB
        # self.ktcBB = 2*self.ktcBB
        # self.ktdBB = 2*self.ktdBB
        
        # Create a new matrix to store the averaged values
        avg_r = r.copy()

        # Calculate the averaged values
        avg_r[1:, 1:] = (r[1:, 1:] + r[1:, :-1]) / 2
        
        # For the first time point, we can just keep the original values as there is no previous point to average with
        avg_r[1:, 0] = r[1:, 0]

        np.savetxt(output_filepath, avg_r, delimiter=',')

    def get_parameters_for_KMC(self, inputs = {}):

        kmc_params = {
            'I_c0': self.I0,
            'A_c0': self.A0,
            'B_c0': self.B0,
            'kd': self.kd,
            'kpAA': self.kpAA,
            'kpAB': self.kpAB,
            'kpBA': self.kpBA,
            'kpBB': self.kpBB,
            'kdAA': self.kdAA,
            'kdAB': self.kdAB,
            'kdBA': self.kdBA,
            'kdBB': self.kdBB,
            'ktcAA': self.ktcAA,
            'ktcAB': self.ktcAB,
            'ktcBB': self.ktcBB,
            'ktdAA': self.ktdAA,
            'ktdAB': self.ktdAB,
            'ktdBB': self.ktdBB
        }
        kmc_params.update(inputs)
        return kmc_params

    def __init__(self, k, f, y0):
        (
            # Initiator decomposition rate constant
            self.kd,

            # Propagation rate constants
            self.kpAA, self.kpAB, self.kpBA, self.kpBB,

            # Depropagation rate constants
            self.kdAA, self.kdAB, self.kdBA, self.kdBB,

            # Termination rate constants
            self.ktcAA, self.ktdAA, 
            self.ktcAB, self.ktdAB, 
            self.ktcBB, self.ktdBB,
        ) = k

        self.kpRA = self.kpAA
        self.kpRB = self.kpBB

        self.ktcAA = 1.0*self.ktcAA
        self.ktdAA = 1.0*self.ktdAA
        self.ktcAB = 1.0*self.ktcAB
        self.ktdAB = 1.0*self.ktdAB
        self.ktcBB = 1.0*self.ktcBB
        self.ktdBB = 1.0*self.ktdBB

        self.ktAA = self.ktcAA + self.ktdAA
        self.ktAB = self.ktcAB + self.ktdAB
        self.ktBB = self.ktcBB + self.ktdBB

        self.k = k
        self.f = f

        # Define initial concentrations
        self.y0 = y0
        self.I0 = self.y0[0]
        self.A0 = self.y0[2]
        self.B0 = self.y0[3]


    def solve(self, t_span, num_points=100, remove_zero=True, solve_chain_model=True,
              solve_sequence_model=True, **kwargs): 
        
        solve_params = {
            'method': 'LSODA',
            'first_step': 1e-12,
            'rtol': 1e-10,
            'atol': 1e-10
        }   
        
        self.solve_chain_model = solve_chain_model
        self.solve_sequence_model = solve_sequence_model
    
        # Solve the method of moments equations
        self.t_span = t_span

        self.mom_sol = solve_ivp(self.kinetic_odes, self.t_span, self.y0, **solve_params)
        
        # Interpolate method of moments solution
        self.t = np.linspace(np.min(self.mom_sol.t), np.max(self.mom_sol.t), num_points)
        self.sol = np.zeros([len(self.mom_sol.y), num_points])

        # print(self.mom_sol.t)
        for i in range(len(self.sol)):
            interp_func = interp1d(self.mom_sol.t, self.mom_sol.y[i], kind='cubic')
            self.sol[i,:] = interp_func(self.t)

        # Remove zero point from self.sol
        if remove_zero:
            self.t = self.t[1:]
            self.sol = self.sol[:,1:]
            
        self.t_s = self.t
        self.t_min = self.t / 60
        self.t_hours = self.t / 3600
            
        self._set_variables(self.sol)
            
    def set_variable_eval(self, t_interp=None, num_points=None):
        
        if num_points is not None:
            t_interp = np.linspace(self.t[0], self.t[-1], num_points)
        if t_interp is None:
            raise ValueError("t_interp must be provided.")
        
        # Confirm that t_interp is within the range of self.t
        if np.min(t_interp) < np.min(self.t) or np.max(t_interp) > np.max(self.t):
            raise ValueError(f"t_interp must be within the range of {np.min(self.t)} to {np.max(self.t)} seconds.")
        
        sol_interp = np.zeros([len(self.sol), len(t_interp)])
        
        for i in range(len(self.sol)):
            interp_func = interp1d(self.t, self.sol[i], kind='cubic')
            sol_interp[i,:] = interp_func(t_interp)
            
        self.sol_ = sol_interp
        self.t_s = t_interp
        self.t_min = t_interp / 60
        self.t_hours = t_interp / 3600
            
        self._set_variables(sol_interp)     

    def _set_variables(self, sol):
        
        # Define monomer concentrations
        eps = 1e-40
        self.cI  = sol[0]
        self.cR  = sol[1]
        self.cA  = sol[2]
        self.cB  = sol[3]
        self.cRA = sol[4] + eps
        self.cRB = sol[5] + eps

        # Define chain 
        self.cPAA = self.PAA_0 = sol[6] + eps
        self.cPAB = self.PAB_0 = sol[7] + eps
        self.cPBA = self.PBA_0 = sol[8] + eps
        self.cPBB = self.PBB_0 = sol[9] + eps
        self.cPD  = self.PD_0 = sol[10] + eps

        self.cPA = self.PA_0 = self.cRA + self.cPAA + self.cPBA
        self.cPB = self.PB_0 = self.cRA + self.cPAB + self.cPBB
        self.cP  = self.aP_0  = self.cPA + self.cPB
        self.P_0 = self.aP_0 + self.PD_0

        # Calculate chian end dyad fractions and set Markov chain end estimations
        self.fPAA = self.pPAAA = self.pPAAB = self.cPAA / self.cPA
        self.fPAB = self.pPABA = self.pPABB = self.cPAB / self.cPB
        self.fPBA = self.pPBAA = self.pPBAB = self.cPBA / self.cPA
        self.fPBB = self.pPBBA = self.pPBBB = self.cPBB / self.cPB

        # Define copolymer composition
        self.CA = (self.A0 - self.cA + eps) / (self.A0 - self.cA + self.B0 - self.cB + eps)
        self.CB = (self.B0 - self.cB + eps) / (self.A0 - self.cA + self.B0 - self.cB + eps)

        # Define chain model moments
        self.PAA_1 = sol[11] + eps
        self.PAB_1 = sol[12] + eps
        self.PBA_1 = sol[13] + eps
        self.PBB_1 = sol[14] + eps
        self.PD_1  = sol[15] + eps

        self.aPA_1 = self.PAA_1 + self.PBA_1
        self.aPB_1 = self.PAB_1 + self.PBB_1
        self.aP_1  = self.aPA_1  + self.aPB_1
        self.P_1   = self.aP_1 + self.PD_1

        self.PAA_2 = sol[16] + eps
        self.PAB_2 = sol[17] + eps
        self.PBA_2 = sol[18] + eps
        self.PBB_2 = sol[19] + eps
        self.PD_2  = sol[20] + eps

        self.PA_2 = self.PAA_2 + self.PBA_2
        self.PB_2 = self.PAB_2 + self.PBB_2
        self.aP_2  = self.PA_2  + self.PB_2
        self.P_2 = self.aP_2 + self.PD_2
        
        # Calculate average chain model parameters
        
        self.nAvgCL = self.NACL = self.P_1 / self.P_0
        self.wAvgCL = self.WACL = self.P_2 / self.P_1
        
        self.PDI  = self.WACL / self.NACL
        self.rpos = self.NACL / max(self.NACL)

        # Define sequence model moments
        self.aSA_0 = sol[21] + eps
        self.aSB_0 = sol[22] + eps
        self.aSA_1 = sol[23] + eps
        self.aSB_1 = sol[24] + eps
        self.aSA_2 = sol[25] + eps
        self.aSB_2 = sol[26] + eps

        self.iSA_0 = sol[27] + eps
        self.iSB_0 = sol[28] + eps
        self.iSA_1 = sol[29] + eps
        self.iSB_1 = sol[30] + eps
        self.iSA_2 = sol[31] + eps
        self.iSB_2 = sol[32] + eps

        self.SA_0 = self.aSA_0 + self.iSA_0
        self.SB_0 = self.aSB_0 + self.iSB_0
        self.SA_1 = self.aSA_1 + self.iSA_1
        self.SB_1 = self.aSB_1 + self.iSB_1
        self.SA_2 = self.aSA_2 + self.iSA_2
        self.SB_2 = self.aSB_2 + self.iSB_2

        # Define average sequence model parameters

        self.aNASL_A = self.aSA_1 / self.aSA_0
        self.aNASL_B = self.aSB_1 / self.aSB_0
        self.aWASL_A = self.aSA_2 / self.aSA_1
        self.aWASL_B = self.aSB_2 / self.aSB_1
        self.aSDI_A  = self.aWASL_A / self.aNASL_A
        self.aSDI_B  = self.aWASL_B / self.aNASL_B
        
        self.aNASL_std_A = np.sqrt((self.aSDI_A - 1) * self.aNASL_A)
        self.aNASL_std_B = np.sqrt((self.aSDI_B - 1) * self.aNASL_B)

        self.iNASL_A = self.iSA_1 / self.iSA_0
        self.iNASL_B = self.iSB_1 / self.iSB_0
        self.iWASL_A = self.iSA_2 / self.iSA_1
        self.iWASL_B = self.iSB_2 / self.iSB_1
        self.iSDI_A  = self.iWASL_A / self.iNASL_A
        self.iSDI_B  = self.iWASL_B / self.iNASL_B
        
        self.iNASL_std_A = np.sqrt((self.iSDI_A - 1) * self.iNASL_A)
        self.iNASL_std_B = np.sqrt((self.iSDI_B - 1) * self.iNASL_B)

        self.nAvgSLA = self.NASL_A = self.SA_1 / self.SA_0
        self.nAvgSLB = self.NASL_B = self.SB_1 / self.SB_0
        self.wAvgSLA = self.WASL_A = self.SA_2 / self.SA_1
        self.wAvgSLB = self.WASL_B = self.SB_2 / self.SB_1
        self.SDI_A  = self.WASL_A / self.NASL_A
        self.SDI_B  = self.WASL_B / self.NASL_B
        
        self.NASL_std_A = np.sqrt((self.SDI_A - 1) * self.NASL_A)
        self.NASL_std_B = np.sqrt((self.SDI_B - 1) * self.NASL_B)


        # Calculate conversion
        self.x_0  = self.P_1 / (self.P_1 + self.cA + self.cB) # Chain model
        self.x_1  = (self.A0 + self.B0 - self.cA - self.cB) / (self.A0 + self.B0) # Concentration
        self.x    = self.x_1
        self.xA   = (self.A0 - self.cA) / self.A0
        self.xB   = (self.B0 - self.cB) / self.B0
        
        self.p_AA = self.kpAA*self.cPA*self.cA
        self.d_AA = self.kdAA*self.cPAA
        self.p_AB = self.kpAB*self.cPA*self.cB
        self.d_AB = self.kdAB*self.cPAB
        self.p_BA = self.kpBA*self.cPB*self.cA
        self.d_BA = self.kdBA*self.cPBA
        self.p_BB = self.kpBB*self.cPB*self.cB
        self.d_BB = self.kdBB*self.cPBB
        
        self.rAA = np.abs(self.p_AA - self.d_AA) + 0*self.fPAA*(self.d_AA + self.d_AB) - 0*self.fPAA*(self.p_AA + self.p_AB)
        self.rAB = np.abs(self.p_AB - self.d_AB) + 0*self.fPAB*(self.d_BB + self.d_BA) - 0*self.fPAB*(self.p_BB + self.p_BA)
        self.rBA = np.abs(self.p_BA - self.d_BA) + 0*self.fPBA*(self.d_AA + self.d_AB) - 0*self.fPBA*(self.p_AA + self.p_AB)
        self.rBB = np.abs(self.p_BB - self.d_BB) + 0*self.fPBB*(self.d_BB + self.d_BA) - 0*self.fPBB*(self.p_BB + self.p_BA)
        
        self.RpA = self.p_AA + self.p_BA - self.d_AA - self.d_BA
        self.RpB = self.p_AB + self.p_BB - self.d_AB - self.d_BB
        

        self.PAA = self.rAA / (self.rAA + self.rBA)
        self.PBA = self.rBA / (self.rAA + self.rBA)
        self.PAB = self.rAB / (self.rAB + self.rBB)
        self.PBB = self.rBB / (self.rAB + self.rBB)
        
        self.pAA = self.rAA / (self.rAA + self.rAB)
        self.pAB = self.rAB / (self.rAA + self.rAB)
        self.pBA = self.rBA / (self.rBA + self.rBB)
        self.pBB = self.rBB / (self.rBA + self.rBB)
        
        self.nAvgAinst = 1/(1-self.pAA)
        self.nAvgBinst = 1/(1-self.pBB)
        
        self.wAvgAinst = (1+self.pAA)/(1-self.pAA)
        self.wAvgBinst = (1+self.pBB)/(1-self.pBB)
        
        def geom_dist(p, n_max):
            return np.array([(1-p)**n * p for n in range(1, n_max+1)])

    # def solve(self, t_span, num_points=20, remove_zero=False):


    #     self.solve_chain_model = True
    #     self.solve_sequence_model = True
        
    #     # Solve the method of moments equations
    #     self.t_span = t_span
    #     # method = None
    #     # if self.kdAA*self.kdAB*self.kdBA*self.kdBB < 1e-5:
    #     #     method = 'BDF'
    #     # else:
    #     #     method = 'LSODA'
    #     # method = 'LSODA'
    #     method = 'LSODA'
    #     # method = 'LSODA'
    #     self.mom_sol = solve_ivp(self.kinetic_odes, self.t_span, self.y0, method=method, rtol=1e-10, atol=1e-10, first_step=1e-10)

    #     # Interpolate method of moments solution
    #     # num_points = 20
    #     self.t = np.linspace(np.min(self.mom_sol.t), np.max(self.mom_sol.t), num_points)
    #     self.t_hours = self.t / 3600
    #     self.sol = np.zeros([len(self.mom_sol.y), num_points])

    #     self.time_points = self.t

    #     for i in range(len(self.sol)):
    #         interp_func = interp1d(self.mom_sol.t, self.mom_sol.y[i], kind='cubic')
    #         self.sol[i,:] = interp_func(self.t)

    #     # Remove zero point from self.sol
    #     if remove_zero:
    #         self.t = self.t[1:]
    #         self.sol = self.sol[:,1:]

    #     # Define monomer concentrations
    #     eps = 1e-40
    #     self.cI = self.sol[0] + eps
    #     self.cR = self.sol[1] + eps
    #     self.cA = self.sol[2] + eps
    #     self.cB = self.sol[3] + eps
    #     self.cRA = self.sol[4] + eps
    #     self.cRB = self.sol[5] + eps

    #     self.cPAA = self.PAA_0 = self.sol[6] + eps
    #     self.cPAB = self.PAB_0 = self.sol[7] + eps
    #     self.cPBA = self.PBA_0 = self.sol[8] + eps
    #     self.cPBB = self.PBB_0 = self.sol[9] + eps
    #     self.cPD  = self.PD_0  = self.sol[10] + eps

    #     self.cPA = self.PA_0 = self.cPAA + self.cPBA
    #     self.cPB = self.PB_0 = self.cPAB + self.cPBB
    #     self.cP = self.P_0 = self.cPA + self.cPB + self.cPD

    #     # Calculate chian end dyad fractions
    #     self.fPAA = self.pPAAA = self.pPAAB = self.cPAA / self.cPA
    #     self.fPAB = self.pPBAA = self.pPBAB = self.cPAB / self.cPB
    #     self.fPBA = self.pPABA = self.pPABB = self.cPBA / self.cPA
    #     self.fPBB = self.pPBBA = self.pPBBB = self.cPBB / self.cPB

    #     # Define copolymer composition
    #     self.CA = (self.A0 - self.cA) / (self.A0 - self.cA + self.B0 - self.cB)
    #     self.CB = (self.B0 - self.cB) / (self.A0 - self.cA + self.B0 - self.cB)

    #     self.x    = (self.A0 + self.B0 - self.cA - self.cB) / (self.A0 + self.B0)
    #     self.xA   = (self.A0 - self.cA) / self.A0
    #     self.xB   = (self.B0 - self.cB) / self.B0
        
    #     self.nAvgCL = None
    #     self.wAvgCL = None
    #     self.nAvgSLA = None
    #     self.wAvgSLA = None
    #     self.nAvgSLB = None
    #     self.wAvgSLB = None

    #     if not self.solve_chain_model:
    #         return

    #     # Define chain model moments
    #     self.PAA_1 = self.sol[11] + eps
    #     self.PAB_1 = self.sol[12] + eps
    #     self.PBA_1 = self.sol[13] + eps
    #     self.PBB_1 = self.sol[14] + eps
    #     self.PD_1  = self.sol[15] + eps

    #     self.PA_1 = self.PAA_1 + self.PBA_1
    #     self.PB_1 = self.PAB_1 + self.PBB_1
    #     self.P_1 = self.PA_1 + self.PB_1 + self.PD_1

    #     self.PAA_2 = self.sol[16] + eps
    #     self.PAB_2 = self.sol[17] + eps
    #     self.PBA_2 = self.sol[18] + eps
    #     self.PBB_2 = self.sol[19] + eps
    #     self.PD_2  = self.sol[20] + eps

    #     self.PA_2 = self.PAA_2 + self.PBA_2
    #     self.PB_2 = self.PAB_2 + self.PBB_2
    #     self.P_2 = self.PA_2 + self.PB_2 + self.PD_2

    #     # Define chain model moments
    #     # self.P_0 = self.sol[6]  + self.sol[7]  + self.sol[8]  + self.sol[9]  + self.sol[10] + 1e-40
    #     # self.P_1 = self.sol[11] + self.sol[12] + self.sol[13] + self.sol[14] + self.sol[15] + 1e-40
    #     # self.P_2 = self.sol[16] + self.sol[17] + self.sol[18] + self.sol[19] + self.sol[20] + 1e-40
    #     # self.PD_0 = self.sol[10] + 1e-40
    #     # self.PD_1 = self.sol[15] + 1e-40
    #     # self.PD_2 = self.sol[20] + 1e-40


    #     # Calculate conversion
    #     self.x_0  = self.P_1 / (self.P_1 + self.cA + self.cB) # Chain model
    #     # self.x_1  = (self.A0 + self.B0 - self.cA - self.cB) / (self.A0 + self.B0) # Concentration
        

    #     # Calculate average chain model parameters
    #     self.NACL = self.P_1 / self.P_0
    #     self.WACL = self.P_2 / self.P_1
    #     self.PDI  = self.WACL / self.NACL

    #     self.nAvgCL = self.NACL
    #     self.wAvgCL = self.WACL

    #     if not self.solve_sequence_model:
    #         return

    #     # Define sequence model moments
    #     self.aSA_0 = self.sol[21] + eps
    #     self.aSB_0 = self.sol[22] + eps
    #     self.aSA_1 = self.sol[23] + eps
    #     self.aSB_1 = self.sol[24] + eps
    #     self.aSA_2 = self.sol[25] + eps
    #     self.aSB_2 = self.sol[26] + eps

    #     self.iSA_0 = self.sol[27] + eps
    #     self.iSB_0 = self.sol[28] + eps
    #     self.iSA_1 = self.sol[29] + eps
    #     self.iSB_1 = self.sol[30] + eps
    #     self.iSA_2 = self.sol[31] + eps
    #     self.iSB_2 = self.sol[32] + eps

    #     self.SA_0 = self.aSA_0 + self.iSA_0
    #     self.SB_0 = self.aSB_0 + self.iSB_0
    #     self.SA_1 = self.aSA_1 + self.iSA_1
    #     self.SB_1 = self.aSB_1 + self.iSB_1
    #     self.SA_2 = self.aSA_2 + self.iSA_2
    #     self.SB_2 = self.aSB_2 + self.iSB_2

    #     # self.nAvgSLA = self.SA_1 / self.SA_0
    #     # self.nAvgSLB = self.SB_1 / self.SB_0

    #     # self.wAvgSLA = self.SA_2 / self.SA_1
    #     # self.wAvgSLB = self.SB_2 / self.SB_1

    #     # self.SDI_A  = self.wAvgSLA / self.nAvgSLA
    #     # self.SDI_B  = self.wAvgSLB / self.nAvgSLB



    #     # self.SA_0 = self.sol[21] + self.sol[27] + 1e-40
    #     # self.SB_0 = self.sol[22] + self.sol[28] + 1e-40
    #     # self.SA_1 = self.sol[23] + self.sol[29] + 1e-40
    #     # self.SB_1 = self.sol[24] + self.sol[30] + 1e-40
    #     # self.SA_2 = self.sol[25] + self.sol[31] + 1e-40
    #     # self.SB_2 = self.sol[26] + self.sol[32] + 1e-40

    #     # Calculate average sequence model parameters
    #     self.NASL_A = self.SA_1 / self.SA_0
    #     self.NASL_B = self.SB_1 / self.SB_0
    #     self.WASL_A = self.SA_2 / self.SA_1
    #     self.WASL_B = self.SB_2 / self.SB_1
    #     self.SDI_A  = self.WASL_A / self.NASL_A
    #     self.SDI_B  = self.WASL_B / self.NASL_B

        


    #     self.nAvgSLA = self.NASL_A
    #     self.nAvgSLB = self.NASL_B
    #     self.wAvgSLA = self.WASL_A
    #     self.wAvgSLB = self.WASL_B

    
    #     # # Solve the Method of Moments system of ODEs
    #     # # method = 'LSODA'
    #     # method = 'BDF'
    #     # # if self.kdAA*self.kdAB*self.kdBA*self.kdBB < 1e-5:
    #     # #     method = 'BDF'
    #     # # else:
    #     # #     method = 'LSODA'

    #     # self.mom_sol = solve_ivp(self.kinetic_odes, self.t_span, self.y0, method=method, rtol=1e-10, atol=1e-10, first_step=1e-10)



    #     # # num_points = 500
    #     # self.t = np.linspace(np.min(self.mom_sol.t), np.max(self.mom_sol.t), num_points)
    #     # self.t_hours = self.t / 3600
    #     # self.sol = np.zeros([len(self.mom_sol.y), num_points])

    #     # for i in range(len(self.sol)):
    #     #     interp_func = interp1d(self.mom_sol.t, self.mom_sol.y[i], kind='cubic')
    #     #     new_y = interp_func(self.t)
    #     #     self.sol[i,:] = new_y

    #     # self.A  = self.sol[2]
    #     # self.B  = self.sol[3]

    #     # self.xA = (self.A0 - self.A) / self.A0
    #     # self.xB = (self.B0 - self.B) / self.B0
    #     # self.x  = (self.A0 + self.B0 - self.A - self.B) / (self.A0 + self.B0)

    #     # # self.save_rxn_rates(self.t, self.sol)

