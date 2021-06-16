'''
Created on Jun 11, 2021

@author: opalazhc
'''
import numpy as np

#constants
Temperature_K = 293 # K

# These are in polynomial for, pulled from a thesis...
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
# Above polynomials are used to evaluate the eq'm constants via polyval  
# from numpy: 
# https://numpy.org/doc/stable/reference/generated/numpy.polyval.html
# Read up on polyval.

#old method, polyval
K_FeOH_hydr = 10 ** (np.polyval(KFeOHPOLYNOMIAL, Temperature_K))
K_FeOH2_hydr = 10 ** (np.polyval(KFeOH2POLYNOMIAL, Temperature_K))
K_FeOH3_hydr = 10 ** (np.polyval(KFeOH3POLYNOMIAL, Temperature_K))

# Evaluating the above constants at a temperature from the Tremaine & LeBlanc
# paper to make sure they are correct 

import math

print('The logarithmic equilibrium constants at a temperature of' , Temperature_K,'K' '\n'
'from the Tremaine & LeBlanc paper using the old method, polyval, for' '\n'
'FeOH, FeOH2 and FeOH3, respectively, are:')
print(math.log10(K_FeOH_hydr), ',' ,math.log10(K_FeOH2_hydr), 'and', 
    math.log10(K_FeOH3_hydr), '.',)

#Making a code using the new method that puts them into a function form
#The input is temperature, the output is the constants

#example of new function, poly1d
#K_FeOH_hydr = np.poly1d(KFeOHPOLYNOMIAL)
#K_FeOH_hydr = 10 ** K_FeOH_hydr(Temperature_K)


while True:
    input_temperature = int(input('\n'
'Input one of the following temperatures (in Kelvin): '
'373, 423, 473, 523, 573' 
'\n'))
    if input_temperature == 373:
        print('The temperature that was inputted was: ', input_temperature,'K')
        break
    elif input_temperature == 423:
        print('The temperature that was inputted was: ', input_temperature,'K')
        break
    elif input_temperature == 473:
        print('The temperature that was inputted was: ', input_temperature,'K')
        break  
    elif input_temperature == 523:
        print('The temperature that was inputted was: ', input_temperature,'K')
        break
    elif input_temperature == 573:
        print('The temperature that was inputted was: ', input_temperature,'K')
        break
    else:
        print('This is an invalid input. Please try again.')

LOG_K_FeOH_hydr = np.poly1d(KFeOHPOLYNOMIAL)
LOG_K_FeOH2_hydr = np.poly1d(KFeOH2POLYNOMIAL)
LOG_K_FeOH3_hydr = np.poly1d(KFeOH3POLYNOMIAL)

K_FeOH_hydr = 10 ** (LOG_K_FeOH_hydr(input_temperature))
K_FeOH_hydr = 10 ** (LOG_K_FeOH_hydr(input_temperature))
K_FeOH_hydr = 10 ** (LOG_K_FeOH_hydr(input_temperature))

print('\n' 'At this temperature, the equilibrium constants' '\n'
    'for FeOH, FeOH2 and FeOH3, respectively,' '\n'
    'from the Tremaine & LeBlanc paper using the new method, poly1d, are:')
    
print(K_FeOH_hydr)
print(K_FeOH2_hydr)
print(K_FeOH3_hydr)