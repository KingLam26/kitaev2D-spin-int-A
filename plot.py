import os
from decimal import Decimal, getcontext
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi

getcontext().prec = 30
spin_S = 2


##### Koga results #####
spin_2_coeff_a = Decimal('18604521')
spin_2_coeff_b = Decimal('-82758048')
spin_2_coeff_c = Decimal('129273554')
spin_2_coeff_d = Decimal('5138022400')

spin_3_coeff_a = Decimal('36052814083126422740')
spin_3_coeff_b = Decimal('-176028114277347622010')
spin_3_coeff_c = Decimal('287126525350219384887')
spin_3_coeff_d = Decimal('152769160756403896320000')

spin_1_5_coeff_a = Decimal('3214648723397092084')
spin_1_5_coeff_b = Decimal('1646995686930432837306')
spin_1_5_coeff_c = Decimal('91522768044989658195')
spin_1_5_coeff_d = Decimal('3076979551468152422400000')

##### extract results #####
"""
def read_results(filename):
    bp_dict = {}
    with open(filename, 'r') as file:
        for line in file:
            # Strip whitespace from the line
            line = line.strip()
            
            # Ignore empty lines
            if not line:
                continue
            
            # Split the line into key and value at the first colon found
            key, value = line.split(':', 1)
            
            # Remove any extra whitespace from key and value
            key = key.strip()
            value = value.strip()
            
            # define the parameters
            if key == "spin_S":
                spin_S = int(value)
            elif key == "multi_factor":
                multi_factor = int(value)
            elif key[:3] == "bp_":
                bp_dict[key] = [int(digit) for digit in str(value)]
    
    return spin_S, multi_factor, bp_dict


spin_S, multi_factor, bp_dict = read_model_params('params-3.txt')



# Print or use the extracted values
for key, value in bp_values.items():
    print(f"{key}: {value}")
"""

##### compute J_eff_sum #####
def trigo(t):
    J = 1
    x = Decimal(J) * Decimal(cos(t))
    y = Decimal(J) * Decimal(sin(t))
    return x,y

def J_eff_sum_lam(t):
    
    x, y = trigo(t)

    # Calculate the sum of terms
    """
    total_sum = Decimal(
        bp_values['bp_1'] * (x**6) * (y**6) +
        bp_values['bp_2'] * (x**6) * (y**6) +
        (bp_values['bp_3'] + bp_values['bp_13']) * (x**8) * (y**4) +
        (bp_values['bp_4'] + bp_values['bp_14']) * (x**4) * (y**8) +
        (bp_values['bp_5'] + bp_values['bp_15']) * (x**8) * (y**4) +
        (bp_values['bp_6'] + bp_values['bp_16']) * (x**4) * (y**8) +
        bp_values['bp_7'] * (x**2) * (y**10) +
        bp_values['bp_17'] * (x**10) * (y**2) +
        (bp_values['bp_8'] + bp_values['bp_18']) * (x**6) * (y**6) +
        (bp_values['bp_9'] + bp_values['bp_19']) * (x**6) * (y**6)
    )

    total_sum = Decimal(
        bp_values['bp_1'] * (x**6) * (y**6) +
        bp_values['bp_2'] * (x**6) * (y**6) +
        (bp_values['bp_3'] + bp_values['bp_13']) * (x**8) * (y**4) +
        (bp_values['bp_4'] + bp_values['bp_14']) * (x**4) * (y**8) +
        (bp_values['bp_5'] + bp_values['bp_15']) * (x**8) * (y**4) +
        (bp_values['bp_6'] + bp_values['bp_16']) * (x**4) * (y**8) +
        bp_values['bp_7'] * (x**2) * (y**10) +
        bp_values['bp_17'] * (x**10) * (y**2) +
        (bp_values['bp_8'] + bp_values['bp_18']) * (x**6) * (y**6) +
        (bp_values['bp_9'] + bp_values['bp_19']) * (x**6) * (y**6)
    )
    """
    
    if spin_S == 2:
        sum = Decimal(
            spin_2_coeff_a * (x**8 + y**8) * 4+
            spin_2_coeff_b * (x**6) * (y**2) * 8 + spin_2_coeff_b * (x**2) * (y**6) * 8+ 
            spin_2_coeff_c * (x**4) * (y**4) * 16
        )
        sum/= spin_2_coeff_d
    
    return sum

def J_eff_sum_koga(t):

    x, y = trigo(t)

    if spin_S == 1.5:
        sum = Decimal(
            spin_1_5_coeff_a * (x**6) * (y**6) +
            spin_1_5_coeff_b * (x**2) * (y**2) * ((x**2) - (y**2))**4 +
            spin_1_5_coeff_c * (x**4) * (y**4) * ((x**2) - (y**2))**2
        )
        sum/=spin_1_5_coeff_d
    
    if spin_S == 2:
    	sum = Decimal(
    	    spin_2_coeff_a * (x**8 + y**8) +
        	spin_2_coeff_b * (x**6) * (y**2) + spin_2_coeff_b * (x**2) * (y**6) + 
        	spin_2_coeff_c * (x**4) * (y**4)
    	)
    	sum/= spin_2_coeff_d
    
    return sum


##### plotting #####
# Generate values of t from 0 to pi/4
t_values = np.linspace(0, pi/4, 5000)
#t_values = np.linspace(0.18, 0.19, 5000) * pi

# Calculate coupling constants for each t
def normalized_values(func, t_values):
    function_values = [Decimal(func(t)) for t in t_values]
    
    # Find maximum magnitude
    max_value = max(abs(val) for val in function_values)
    
    # Normalize the function values
    normalized_values = [val / max_value for val in function_values]
    
    return normalized_values

normalized_values_lam = normalized_values(J_eff_sum_lam, t_values)
normalized_values_koga = normalized_values(J_eff_sum_koga, t_values)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t_values / pi, normalized_values_lam, label='lam', linewidth = 2)
plt.plot(t_values / pi, normalized_values_koga, label='koga', linewidth = 1)
plt.xlabel('t (radians) / pi')
plt.ylabel('Nomalized J_eff_sum')
plt.title('Nomalized J_eff_sum against theta')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.savefig('plot.svg', format='svg')