import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

class PenultimateModel():

    def kinetic_odes(self, t, y):
        dydt = np.zeros_like(y)

        # Define small molecule concentrations
        I  = y[0]
        R  = y[1]
        A  = y[2]
        B  = y[3]
        RA = y[4]
        RB = y[5]

        # Define chain end dyad fractions
        PAA_0 = y[6]
        PAB_0 = y[7]
        PBA_0 = y[8]
        PBB_0 = y[9]

        PA_0 = PAA_0 + PBA_0 + 1e-40
        PB_0 = PAB_0 + PBB_0 + 1e-40

        # Calculate chain end dyad fractions
        fPAA = PAA_0 / PA_0
        fPAB = PAB_0 / PB_0
        fPBA = PBA_0 / PA_0
        fPBB = PBB_0 / PB_0

        # Define chain end triad fractions from Markov Chain factorization
        pPAAA = pPAAB = fPAA
        pPBAA = pPBAB = fPBA
        pPABA = pPABB = fPAB
        pPBBA = pPBBB = fPBB


        #####################################
        ### Small molecule concentrations ###
        #####################################

        # Initiator concentration
        dydt[0] = -self.kd*I

        # Free radical concentration
        dydt[1] = 2*self.f*self.kd*I - self.kpRA*R*A - self.kpRB*R*B

        # A monomer concentration
        dydt[2] = -A*(self.kpRA*R + self.kpAA*(RA+PA_0) + self.kpBA*(RB+PB_0)) + self.kdAA*PAA_0 + self.kdBA*PBA_0

        # B monomer concentration
        dydt[3] = -B*(self.kpRB*R + self.kpBB*(RB+PB_0) + self.kpAB*(RA+PA_0)) + self.kdBB*PBB_0 + self.kdAB*PAB_0

        # RA dyad concentration
        dydt[4] = self.kpRA*R*A - self.kpAA*RA*A - self.kpAB*RA*B

        # RB dyad concentration
        dydt[5] = self.kpRB*R*B - self.kpBB*RB*B - self.kpBA*RB*A

        ########################################
        ### Chain model (0th order moments) ####
        ########################################

        # PAA concentration
        dydt[6]   = self.kpAA*RA*A + self.kpAA*PA_0*A - \
                    self.kpAA*PAA_0*A - self.kpAB*PAA_0*B - \
                    self.kdAA*PAA_0 + self.kdAA*pPAAA*PAA_0 + self.kdAB*pPAAB*PAB_0 - \
                    self.ktAA*PAA_0*PA_0 - self.ktAB*PAA_0*PB_0
        # PAB concentration
        dydt[7]   = self.kpAB*RA*B + self.kpAB*PA_0*B - \
                    self.kpBA*PAB_0*A - self.kpBB*PAB_0*B - \
                    self.kdAB*PAB_0 + self.kdBA*pPABA*PBA_0 + self.kdBB*pPABB*PBB_0 - \
                    self.ktAB*PAB_0*PA_0 - self.ktBB*PAB_0*PB_0

        # PBA concentration
        dydt[8]   = self.kpBA*RB*A + self.kpBA*PB_0*A - \
                    self.kpAB*PBA_0*B - self.kpAA*PBA_0*A - \
                    self.kdBA*PBA_0 + self.kdAB*pPBAB*PAB_0 + self.kdAA*pPBAA*PAA_0 - \
                    self.ktAB*PBA_0*PB_0 - self.ktAA*PBA_0*PA_0

        # PBB concentration
        dydt[9]   = self.kpBB*RB*B + self.kpBB*PB_0*B - \
                    self.kpBB*PBB_0*B - self.kpBA*PBB_0*A - \
                    self.kdBB*PBB_0 + self.kdBB*pPBBB*PBB_0 + self.kdBA*pPBBA*PBA_0 - \
                    self.ktBB*PBB_0*PB_0 - self.ktAB*PBB_0*PA_0

        # D concentration
        dydt[10]  = 0.5*self.ktcAA*PA_0*PA_0 + 0.5*self.ktcAB*PA_0*PB_0 + 0.5*self.ktcBB*PB_0*PB_0 + \
                    self.ktdAA*PA_0*PA_0 + 0.5*self.ktdAB*(PA_0*PB_0 + PB_0*PA_0) + self.ktdBB*PB_0*PB_0
        
        ########################################
        ### Chain model (1st order moments) ####
        ########################################

        ########################################
        ### Chain model (2nd order moments) ####
        ########################################

        ###################################################
        ### Sequence model, active (0th order moments) ####
        ###################################################

        ###################################################
        ### Sequence model, active (1st order moments) ####
        ###################################################

        #####################################################
        ### Sequence model, inactive (0th order moments) ####
        #####################################################

        #####################################################
        ### Sequence model, inactive (1st order moments) ####
        #####################################################


        return dydt

    def get_rxn_strings(self):
        # Define all polymer species
        polymer_species = ["PR", "PRA", "PRB", "PAA", "PAB", "PBA", "PBB", "PD"]
        
        # 35 reactions
        reaction_strings = [
            "AIBN -kd-> PR",

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

            "PAA -kdAA-> PAA",
            "PAB -kdAB-> PAB",
            "PBA -kdBA-> PBA",
            "PBB -kdBB-> PBB",

            "PAA -ktcAA-> PD",
            "PBA -ktcAA-> PD",

            "PAA -ktcAB-> PD",
            "PBA -ktcAB-> PD",
            "PAB -ktcAB-> PD",
            "PBB -ktcAB-> PD",

            "PAB -ktcBB-> PD",
            "PBB -ktcBB-> PD",

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

        # self.ktcAA = 0.5*self.ktcAA
        # self.ktdAA = 0.5*self.ktdAA
        # self.ktcAB = 0.5*self.ktcAB
        # self.ktdAB = 0.5*self.ktdAB
        # self.ktcBB = 0.5*self.ktcBB
        # self.ktdBB = 0.5*self.ktdBB

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

        # Terminati7000on by comb.: PA + PAB -> D
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

        # Termination by disp.: PAA + PA -> 2D
        r[28,:] = self.ktdAA*PAA*PA/1.0

        # Termination by disp.: PBA + PA -> 2D
        r[29,:] = self.ktdAA*PBA*PA/1.0

        # Termination by disp.: PAA + PB -> 2D
        r[30,:] = self.ktdAB*PAA*PB/1.0

        # Termination by disp.: PBA + PB -> 2D
        r[31,:] = self.ktdBB*PAB*PB/1.0

        # Termination by disp.: PBB + PB -> 2D
        r[35,:] = self.ktdBB*PBB*PB/1.0

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

        np.savetxt(output_filepath, r, delimiter=',')





    def __init__(self, k, f, y0, t_span):
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

        self.ktcAA = 2*self.ktcAA
        self.ktdAA = 2*self.ktdAA
        self.ktcAB = 2*self.ktcAB
        self.ktdAB = 2*self.ktdAB
        self.ktcBB = 2*self.ktcBB
        self.ktdBB = 2*self.ktdBB

        self.ktAA = self.ktcAA + self.ktdAA
        self.ktAB = self.ktcAB + self.ktdAB
        self.ktBB = self.ktcBB + self.ktdBB

        self.k = k
        self.f = f
        self.y0 = y0
        self.t_span = t_span
        self.I0 = self.y0[0]
        self.A0 = self.y0[2]
        self.B0 = self.y0[3]

        # Solve the Method of Moments system of ODEs
        method = 'LSODA'
        # if self.kdAA*self.kdAB*self.kdBA*self.kdBB < 1e-5:
        #     method = 'BDF'
        # else:
        #     method = 'LSODA'
        self.mom_sol = solve_ivp(self.kinetic_odes, self.t_span, self.y0, method=method, rtol=1e-10, atol=1e-10, first_step=1e-10)


        num_points = 20
        self.t = np.linspace(np.min(self.mom_sol.t), np.max(self.mom_sol.t), num_points)
        self.t_hours = self.t / 3600
        self.sol = np.zeros([len(self.mom_sol.y), num_points])

        for i in range(len(self.sol)):
            interp_func = interp1d(self.mom_sol.t, self.mom_sol.y[i], bounds_error=False, kind='cubic')
            new_y = interp_func(self.t)
            self.sol[i,:] = new_y

        self.A  = self.sol[2]
        self.B  = self.sol[3]

        self.xA = (self.A0 - self.A) / self.A0
        self.xB = (self.B0 - self.B) / self.B0
        self.x  = (self.A0 + self.B0 - self.A - self.B) / (self.A0 + self.B0)

        # self.save_rxn_rates(self.t, self.sol)

