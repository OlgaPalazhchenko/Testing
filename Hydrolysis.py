'''
Created on Jun 11, 2021

@author: opalazhc
'''
import numpy as np


#constants
Temperature_K = 293 # K


# these solve for the log form of the hydrolysis constants for iron,
# e.g.,

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
    
KFeOHPOLYNOMIAL = [-4E-10, 8.1013E-07, -0.00062963, 0.230763, -41.5545]
KFeOH2POLYNOMIAL = [-4E-10, 8.63467E-07, -0.00073731, 0.311499, -67.8248]
KFeOH3POLYNOMIAL = [
    -4.667E-10, 1.0496E-06, -0.000935775, 0.413186, -97.4709
    ]


' 1. add the hydrolysis expressions into the function skeleton provided below' 
' output (returns) should be the hydrolysis constants'
' 2. call the function at a few different temperatures' 
' (no need to write fancy print output)'
'3. make a list/array of temperatures to test and iterate calling function through it' 
'4. make an empty array to save (append) the eventual solutions of the function'
'5. loop through the input temperature list and append solutions to just one of' 
' the hyrolysis constants to it'
'6. print the list' 

#Steps 1 and 2

Temp = 293

def iron_hydrolysis_constants(Temp):
    
    Exponents_K_FeOH_hydr = np.poly1d(KFeOHPOLYNOMIAL)
    Exponents_K_FeOH2_hydr = np.poly1d(KFeOH2POLYNOMIAL)
    Exponents_K_FeOH3_hydr = np.poly1d(KFeOH3POLYNOMIAL)
    
    Log_K_FeOH_hydr = Exponents_K_FeOH_hydr(Temp)
    Log_K_FeOH2_hydr = Exponents_K_FeOH2_hydr(Temp)
    Log_K_FeOH3_hydr = Exponents_K_FeOH3_hydr(Temp)

    K_FeOH_hydr = 10 ** Log_K_FeOH_hydr
    K_FeOH2_hydr = 10 ** Log_K_FeOH2_hydr
    K_FeOH3_hydr = 10 ** Log_K_FeOH3_hydr

    return K_FeOH_hydr, K_FeOH2_hydr, K_FeOH3_hydr

print(iron_hydrolysis_constants(Temp))




#Step 3

TempArray =np.array([293, 323, 373, 423, 473, 523])

def iron_hydrolysis_constants(TempArray):
    
  # for iron_hydrolysis_constants in Temp:
  #     print(iron_hydrolysis_constants)
    
    Exponents_K_FeOH_hydr = np.poly1d(KFeOHPOLYNOMIAL)
    Exponents_K_FeOH2_hydr = np.poly1d(KFeOH2POLYNOMIAL)
    Exponents_K_FeOH3_hydr = np.poly1d(KFeOH3POLYNOMIAL)
    
    Log_K_FeOH_hydr = Exponents_K_FeOH_hydr(TempArray)
    Log_K_FeOH2_hydr = Exponents_K_FeOH2_hydr(TempArray)
    Log_K_FeOH3_hydr = Exponents_K_FeOH3_hydr(TempArray)

    K_FeOH_hydr = 10 ** Log_K_FeOH_hydr
    K_FeOH2_hydr = 10 ** Log_K_FeOH2_hydr
    K_FeOH3_hydr = 10 ** Log_K_FeOH3_hydr

    constants_array = np.array([K_FeOH_hydr,  K_FeOH2_hydr, K_FeOH_hydr])

    Value = np.dot(constants_array, TempArray)

    return Value

#[K_FeOH_hydr, K_FeOH2_hydr, K_FeOH3_hydr]

print(iron_hydrolysis_constants(TempArray))