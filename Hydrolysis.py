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

# Evaluate the above constants at a temperature from the Tremaine & LeBlanc
# paper to make sure they are correct 
# (they are...but you should be able to show this)

# play around with more recent poly1d method...you can find documentation
# for it using above link, just search for poly1d

# example of new function, poly1d
K_FeOH_hydr = np.poly1d(KFeOHPOLYNOMIAL)
K_FeOH_hydr = 10 ** K_FeOH_hydr(Temperature_K)

# Code up the constats using new method and put them all into a function form
# input sould be temperature, output should be the constants

