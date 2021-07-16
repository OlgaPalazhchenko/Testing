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

Temperatures_Kelvin = [293, 323, 373, 423, 473, 523]


def iron_hydrolysis_constants(Temperature):
    
    ' input here is ANY temperature, not hard-coded for just that one array'
    
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


print (iron_hydrolysis_constants(Temperatures_Kelvin)[0], '\n'\
       'array input direct into func')


' looping through :'
' I actually didnt know that python functions loop over an input list'
'and use each element of the list as a separate argument to the function'

' I worry that could get confusing if there are multiple inputs and they are all lists...'
'This is another way to do it:'


' empty array' 
K_Fe2_at_Ts = []

for Temp in Temperatures_Kelvin:
    K_Fe2 = iron_hydrolysis_constants(Temp)[0]
    K_Fe2_at_Ts.append(K_Fe2)

print ('\n',
    K_Fe2_at_Ts, '\n loopig through each item in the input array')




a = [1,2,3]
# b = [0,1,4]


' This is what I mean by having multiple lists as input...there are list on list \
operations (multiplying, adding, etc.) that will not work unless specifically \
coded to do so... but that is a cool thing to know for simple functions with list \
input! '

def test(x):
    answer0 = 2 * a
    
    # this operation doesn't work by multiplying 2 by each list item
    # instead it duplicates list 
    # list input worked for your original function becaus of the np.polyval
    # function.... numpy arrays are designed for those calcs
    # I would recommend doing the looping through and appending method instead
    print (answer0)
    return answer0      

test (a)







