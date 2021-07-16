'''
Created on Jul 16, 2021

@author: opalazhc
'''

import numpy as np
import sympy as sym


SULPHTOTAL = 0.008 # mol/L
FeII_TOTAL = 0.008 # mol/L
FeIII_TOTAL = 0.00001 # mol/L

Temperature_K = 293 # K

derivative = []
ImpurityType = "iron"

# SO42- + H2O <-> HSO4- + OH-
Kb_HSO4 = 10 ** (-12)
K_FeOH_2plus_hydr = 10 ** (-2.19)
K_FeOHplus_hydr = 10 ** (-5.76)
K_FeOH3_0_hydr = 10 ** (-14.30)
K_FeOH4_minus_hydr = 10 ** (-21.79)


def iron_hydrolysis_constants(Temperature):
    
    KFeOHPOLYNOMIAL = [-4E-10, 8.1013E-07, -0.00062963, 0.230763, -41.5545]
    KFeOH2POLYNOMIAL = [-4E-10, 8.63467E-07, -0.00073731, 0.311499, -67.8248]
    KFeOH3POLYNOMIAL = [
    -4.667E-10, 1.0496E-06, -0.000935775, 0.413186, -97.4709
    ]
    
    # log_10 K(T) polynomial functions
    K_FeOH_hydr_plynom = np.poly1d(KFeOHPOLYNOMIAL)
    K_FeOH2_hydr_polynom = np.poly1d(KFeOH2POLYNOMIAL)
    K_FeOH3_hydr_polynom = np.poly1d(KFeOH3POLYNOMIAL)
    
    ' I do not love the Exponents verbage out front...that does not really mean anything'
    'to me ... changed to polynom, as these are polynomials'
    
    Log_K_FeOH_hydr = K_FeOH_hydr_plynom(Temperature)
    Log_K_FeOH2_hydr = K_FeOH2_hydr_polynom(Temperature)
    Log_K_FeOH3_hydr = K_FeOH3_hydr_polynom(Temperature)
    
    # Fe2+ hydrolysis constants as function of input temp
    K_FeOH_hydr = 10 ** Log_K_FeOH_hydr
    K_FeOH2_hydr = 10 ** Log_K_FeOH2_hydr
    K_FeOH3_hydr = 10 ** Log_K_FeOH3_hydr
    

    return K_FeOH_hydr, K_FeOH2_hydr, K_FeOH3_hydr


K_FeOH_hydr, K_FeOH2_hydr, K_FeOH3_hydr = iron_hydrolysis_constants(Temperature_K)
 


def KW_DBHUCK(Temp_K):
    
    # Dissociation constant of water as a function of temperature
    T_Celsius = Temp_K - 298.15
    
    KwPOLYNOMIAL = [
        -7.226E-10, 1.32661E-06, -0.000959311, 0.32765297, -55.86334915
        ]
    KW = 10 ** (np.polyval(KwPOLYNOMIAL, Temp_K))
    
    DEBYE_HUCKEL_POLYNOMIAL = [3.29E-10, -1.709E-07, 0.00003315, -0.0009028, 0.5027]

    DebyeHuckel = (np.polyval(DEBYE_HUCKEL_POLYNOMIAL, T_Celsius))

    return KW, DebyeHuckel

KW, DebyeHuckelConstant = KW_DBHUCK(Temperature_K)


def Ionic_Strength(InputIons, InputCharges):
    # used with hydrolysis function to iterate activity coeffs / hydrolysis distrib
    
    IS = []
    for x, y in zip(InputCharges, InputIons):
        z = (x **2) * y
        IS.append(z)
    # sum of concentrations and charges squared     
    IonicStr = sum(IS)
    
    # choice of ionic strength equation for activity coefficient calculation
    # based on magnitude of ionic strength
    if IonicStr < 0.01:
        IS_eq = ((IonicStr ** 0.5) / (1 + IonicStr ** 0.5)) - 0.2 * IonicStr
    else:
        
        B = 0.3273 
        rad_ion = 4.5e-8  # ionic radius, estimated as that for Fe2+
        
        IS_eq = ((IonicStr ** 0.5) / (1 + B * rad_ion * IonicStr ** 0.5))               
   
    # depends on ionic charge (coeffs 1-3 for absolute charges of 1-3):
    ActCoeff1 = 10 ** (-DebyeHuckelConstant * (1 ** 2) * IS_eq)
    ActCoeff2 = 10 ** (-DebyeHuckelConstant * ((-2) ** 2) * IS_eq)
    ActCoeff3 = 10 ** (-DebyeHuckelConstant * ((-3) ** 2) * IS_eq)
    
    return ActCoeff1, ActCoeff2, ActCoeff3


def pH_hydrolysis_calculator_FeSystem(
        pHinput, i, Conc_H, ActCoeff1, ActCoeff2, ActCoeff3,
        Solution_System):
    # pH calculator mode:
    # at i=0, derivative eqn'n in symbol form, no numerical sol'n out
    # i >0, hydrolysis distrib. recalculated based on H+ input and activity coeff's 
    # Activity coeffs recalculated in IS function based on new H+ 
    # and hydrolysis []'s
    # New H+ found, iterated until convergence ==> pH solver
    
    # hydrolysis distribution mode:
    # Use known H+ to find hydtolysis distribution
    # Activity coeffs iteratatively recalculated in IS function based on known H+
    # ==> iron species concentration and activity coeffs solver
    
    #iteration loops for either of above methods are external to this function
    
    if pHinput == 'yes' and i == 0:
        i = 1
        print (
            'i = 0 is for algebraic solution to H+. \
            Known pH has been selected and i has been changed to 1'
            )
        
    # first iteration, i=0, only establishes derivative expression in symbol 
    # form, no numerical solution
    if i == 0:
        H = sym.Symbol('H')
    else:
        H = Conc_H
    
    if Solution_System == 'FeSO4':
        TotalSulphate = SULPHTOTAL
        SO4_2 = (
            TotalSulphate / (
                1
                + (Kb_HSO4 * ActCoeff2 * H * ActCoeff1 / ActCoeff1)
                )
            )

        HSO4_minus = Kb_HSO4 * SO4_2 * ActCoeff2 * H # * ActCoeff1 / ActCoeff1
        Li = 0
        OH = KW / (H * (ActCoeff1 ** 2))
    
    elif Solution_System == 'LiOH':
        Li = KW / Conc_H # * ActCoeff1 / ActCoeff1
        OH = KW / Conc_H # * ActCoeff1 / ActCoeff1
        SO4_2 = 0
        HSO4_minus = 0
        
    else:
        None
        
    # Fe2 + H2O = Fe(OH)+ + H+ 
    # K = Fe(OH)*a1 + * H* a1 / Fe2*a2
    # FeOH*a1 = K * Fe*a2 / H*a1
    # Fe(OH)*a1 / Fe2+ * a2  = K / H * a1
    
    # Fe2+ + 2H2O = Fe(OH)2 0 + 2H+
    # K = Fe(OH)2 0 (H* a1) ** 2 / Fe2+ * a2
    # Fe(OH)2 0 = K Fe2+ * a2 / (H * a1) ** 2
    # Fe(OH)2 0/ Fe2+ *a2 = K/ (H*a1) ** 2
    
    # Fe2+ + 3H2O = Fe(OH)3 - + 3H+
    # K = Fe(OH)3 - (H* a1) ** 3 / Fe2+ * a2
    # Fe(OH)3 - = K Fe2+ * a2 / (H * a1) ** 3
    # Fe(OH)3 -/ Fe2+ *a2 = K/ (H*a1) ** 3
    
    #FeTOT = Fe2 + FeOH+ + FeOH2 + FeOH3-
    # Fe2 = FeTOT - FeOH+ - FeOH2 - FeOH3-
    # 1+ FeOH+/Fe2 - Fe(OH)2/Fe2 - Fe(OH)3-/Fe2  = FeTot / Fe2
    
    Total_Fe_II = FeII_TOTAL
    Total_Fe_III = FeIII_TOTAL
    
    Fe2 = (
        Total_Fe_II / (
            1
            + (K_FeOH_hydr * ActCoeff2 / (H * ActCoeff1 ** 2))
            + (K_FeOH2_hydr * ActCoeff2 / ((H ** 2) * ActCoeff1 ** 2))
            + (K_FeOH3_hydr * ActCoeff2 / ((H ** 3) * ActCoeff1 ** 4))
            )
        )
    
    FeOH = (
        K_FeOH_hydr * ActCoeff2 * Fe2 / (H * ActCoeff1 ** 2)
        ) 
    FeOH2_0 = (
        K_FeOH2_hydr * ActCoeff2 * Fe2 / ((H ** 2) * ActCoeff1 ** 2)
        )     
    FeOH3_minus = (
        K_FeOH3_hydr * ActCoeff2 * Fe2 / ((H ** 3) * ActCoeff1 ** 4)
        )
    
    # Fe3+ + H2O = Fe(OH)^ 2+ + H+ 
    # K = [Fe(OH)^ 2+]A2 * [H+]A1 / [Fe3+]A3
    # Fe(OH)^ 2+/Fe3 = K * A3 / H+ * A1 * A2
    
    # Fe3+ + 2H2O = Fe(OH)2 ^ 1 + 2H+ 
    # K = [(FeOH)2 +]A2 * [H+]A1 / [Fe3+]A3
    # Fe(OH)2 +/Fe3 = K * A3 / (H+ A1) **2 * A1
    
    # Fe3+ + 3H2O = Fe(OH)30 + 3H+ 
    # K = [FeOH30]A1 * ([H+]A1) **3 / [Fe3+]A3
    # FeOH30/Fe3 A3 = K / (H+ *A1) **3
    
    # Fe3+ + 4H2O = Fe(OH)4 - + 4H+ 
    # K = [FeOH4-]A1 * ([H+]A1)**4 / [Fe3+]A3
    # Fe(OH)4-*a1/Fe3*a3 = K  / (H+)**4 A1**4
    
    Fe3 = (
    Total_Fe_III / (
        1
        + (K_FeOH_2plus_hydr * ActCoeff3 / (H * ActCoeff1 * ActCoeff2))
        + (K_FeOHplus_hydr * ActCoeff3  / ((H ** 2) * ActCoeff1 ** 3))
        + (K_FeOH3_0_hydr * ActCoeff3  / ((H ** 3) * ActCoeff1 ** 3))
        + (K_FeOH4_minus_hydr * ActCoeff3  / ((H ** 4) * ActCoeff1 ** 5))
        )
    )
    
    FeOH_2plus = (
        K_FeOH_2plus_hydr * ActCoeff3 * Fe3 / (H * ActCoeff1 * ActCoeff2)
        )
    FeOH2_plus = (
        K_FeOHplus_hydr * ActCoeff3 * Fe3 / ((H ** 2) * ActCoeff1 ** 3)
        )     
    FeOH3_0 = (
        K_FeOH3_0_hydr * ActCoeff3 * Fe3 / ((H ** 3) * ActCoeff1 ** 3)
        )
    FeOH4_minus = (
        K_FeOH4_minus_hydr * ActCoeff3 * Fe3 / ((H ** 4) * ActCoeff1 ** 5)
        )
    
    # CB: 
    # [H+] + [FeOH+] + 2[Fe2+] + 3[Fe3+] + 2[Fe(OH)^2+] + [Fe(OH)2^+] 
    # = 2[SO42-] + [OH-] + [FeOH3-] + [HSO4-]
    
    
    AllAnions = OH + FeOH3_minus + FeOH4_minus + SO4_2 + HSO4_minus
    
    AllCations = H + 2 * Fe2 + 3 * Fe3 + FeOH + FeOH2_plus + 2 * FeOH_2plus \
    + Li
    
    InputIons = [
        H, OH, Fe2, FeOH, FeOH3_minus, Fe3, FeOH_2plus, FeOH2_plus,
        FeOH4_minus, SO4_2, HSO4_minus, Li
        ]
    
    InputCharges = [1, -1, 2, 1, -1, 3, 2, 1, -1, -2, 1]

    
    FM =  AllAnions - AllCations
    
    if pHinput == 'yes':
        
        print(
            'Fe2:', 100 * Fe2/(Total_Fe_II + Total_Fe_III), 
            'FeOH:',  100 * FeOH / (Total_Fe_II + Total_Fe_III), 
            'FeOH20',  100 * FeOH2_0 / (Total_Fe_II + Total_Fe_III), 
            'FeOH3_minus: ', 100 * FeOH3_minus/ (Total_Fe_II + Total_Fe_III)
            )
        
        return InputIons, InputCharges, Fe2, FeOH, FeOH3_minus, SO4_2, Fe3,\
            FeOH_2plus, FeOH2_plus, FeOH2_0, Li
       
    
    else:
        # sets up equation for H+ as function of all eq'm iron species variables
        if i == 0:    
            DM = FM.diff(H)
            derivative.append(DM)
            RE = 1e-7
            NM = Conc_H
        else:
            
            DM_eval = derivative[0].evalf(subs={'H': Conc_H})
            
            NM = H - FM / DM_eval
            RE = abs((NM - H) / NM)
    
    
  
    # print (
    #     'Fe3:', Fe3/(Total_Fe_II + Total_Fe_III), 
    #     'Fe(OH)^2+:', FeOH_2plus / (Total_Fe_II + Total_Fe_III), 
    #     'Fe(OH2)^+', FeOH2_plus / (Total_Fe_II + Total_Fe_III), 
    #     'FeOH3_0:', FeOH3_0 / (Total_Fe_II + Total_Fe_III),
    #     'FeOH4_-:',  FeOH4_minus/ (Total_Fe_II + Total_Fe_III),
    #
    #     )

        return (
            InputIons, InputCharges, Fe2, FeOH, FeOH3_minus, SO4_2, Fe3,
            FeOH_2plus, FeOH2_plus, FeOH2_0, RE, NM
            )


def hydrolysis_distrib(pH_value, Solution_System):
    
    #desired pH
    a_H = 10 ** (-pH_value)
    
    # initial guesses
    ActCoeff1 = 1
    ActCoeff2 = 1
    ActCoeff3 = 1
   
    # iteratively solve activity coeffs and hydrolysis distrib 
    for i in range(1, 50):
       
        A1 = ActCoeff1
        A2 = ActCoeff2
        A3 = ActCoeff3
        
        
        InputIons, InputCharges, Fe2, FeOH, FeOH3_minus, SO4_2, Fe3,\
            FeOH_2plus, FeOH2_plus, FeOH2_0, Li = \
            pH_hydrolysis_calculator_FeSystem(
            'yes', i, a_H, A1, A2, A3, Solution_System
            )
        
        ActCoeff1, ActCoeff2, ActCoeff3 = Ionic_Strength(
                InputIons, InputCharges
                )
        RE = abs((A1 - ActCoeff1) / A1)
        
        if RE < 1e-7:
            iterations = i
            print (iterations)
            break
            
    
    return Fe2, FeOH, FeOH3_minus, Fe3, FeOH_2plus, FeOH2_plus, FeOH2_0, SO4_2

    
hydrolysis_distrib(pH_value = 10.4, Solution_System = 'LiOH')

' output here is percentage'
