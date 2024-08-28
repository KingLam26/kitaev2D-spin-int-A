import os
from decimal import Decimal, getcontext
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi

getcontext().prec = 30
spin_S = 2
consider_all_combis_switch = True


##### Koga results #####
spin_2_coeff_mu = Decimal('18604521')
spin_2_coeff_nu = Decimal('-82758048')
spin_2_coeff_la = Decimal('129273554')
spin_2_coeff_d = Decimal('5138022400')

spin_3_coeff_mu = Decimal('36052814083126422740')
spin_3_coeff_nu = Decimal('-176028114277347622010')
spin_3_coeff_la = Decimal('287126525350219384887')
spin_3_coeff_d = Decimal('152769160756403896320000')

spin_1_5_coeff_a = Decimal('3214648723397092084')
spin_1_5_coeff_b = Decimal('1646995686930432837306')
spin_1_5_coeff_c = Decimal('91522768044989658195')
spin_1_5_coeff_d = Decimal('3076979551468152422400000')

##### Lam results #####
spin_2_coeff_a = Decimal('18604521')
spin_2_coeff_b = Decimal('92064512')
spin_2_coeff_c = Decimal('-41379024')

spin_3_coeff_a = Decimal('36052814083126422740')
spin_3_coeff_b = Decimal('427101825544440584157')
spin_3_coeff_c = Decimal('124066871221800233745')

spin_4_coeff_a = 
spin_4_coeff_b = 
spin_4_coeff_c = 
spin_4_coeff_d = 
spin_4_coeff_e = 
spin_4_coeff_f = 

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
    if spin_S == 2:
        if consider_all_combis_switch == True:
            x8 = Decimal(spin_2_coeff_a * 4 * (x**8))
            y8 = Decimal(spin_2_coeff_a * 4 * (y**8))
            
            x4_y4 = Decimal((2*spin_2_coeff_a * 4 + spin_2_coeff_b * 16) * ((x**4) * (y**4)))
            
            x6_y2 = Decimal(2*spin_2_coeff_c * 8 * ((x**6) * (y**2)))
            x2_y6 = Decimal(2*spin_2_coeff_c * 8 * ((x**2) * (y**6)))
            
        elif consider_all_combis_switch == False:
            x8 = Decimal(spin_2_coeff_a * (x**8))
            y8 = Decimal(spin_2_coeff_a * (y**8))
            
            x4_y4 = Decimal((2*spin_2_coeff_a + spin_2_coeff_b) * ((x**4) * (y**4)))
            
            x6_y2 = Decimal(2*spin_2_coeff_c * ((x**6) * (y**2)))
            x2_y6 = Decimal(2*spin_2_coeff_c * ((x**2) * (y**6)))

        sum = Decimal(x8 + y8 + x4_y4 + x6_y2 + x2_y6)
        sum/= spin_2_coeff_d

    if spin_S == 3:
        if consider_all_combis_switch == True:
            x12 = Decimal(spin_3_coeff_a * (x**12) * 4)
            y12 = Decimal(spin_3_coeff_a * (y**12) * 4)

            x6_y6 = Decimal((-2*spin_3_coeff_a * 4 - 2*spin_3_coeff_b * 16) * ((x**6) * (y**6)))

            x4_y8 = Decimal((2*spin_3_coeff_c * 8 + spin_3_coeff_b * 16) * (x**4 * y**8))
            x8_y4 = Decimal((2*spin_3_coeff_c * 8 + spin_3_coeff_b * 16) * (y**4 * x**8))

            x2_y10 = Decimal((-2*spin_3_coeff_c) * 8 * (y**10 * x**2))
            x10_y2 = Decimal((-2*spin_3_coeff_c) * 8 * (y**2 * x**10))
        
        elif consider_all_combis_switch == False:
            x12 = Decimal(spin_3_coeff_a * (x**12))
            y12 = Decimal(spin_3_coeff_a * (y**12))
            x6_y6 = Decimal((-2*spin_3_coeff_a - 2*spin_3_coeff_b) * ((x**6) * (y**6)))
            x4_y8 = Decimal((2*spin_3_coeff_c + spin_3_coeff_b) * (x**4 * y**8))
            x8_y4 = Decimal((2*spin_3_coeff_c + spin_3_coeff_b) * (y**4 * x**8))
            x2_y10 = Decimal((-2*spin_3_coeff_c) * (y**10 * x**2))
            x10_y2 = Decimal((-2*spin_3_coeff_c) * (y**2 * x**10))

        sum = Decimal(x12 + y12 + x6_y6 + x4_y8 + x8_y4 + x2_y10 + x10_y2)
        sum = Decimal(sum / spin_3_coeff_d)
    return sum

def J_eff_sum_koga(t):

    x, y = trigo(t)

    if spin_S == 1.5:
        sum = Decimal(spin_1_5_coeff_a * (x**6) * (y**6) +
                      spin_1_5_coeff_b * (x**2) * (y**2) * ((x**2) - (y**2))**4 +
                      spin_1_5_coeff_c * (x**4) * (y**4) * ((x**2) - (y**2))**2)
        sum/=spin_1_5_coeff_d
    
    if spin_S == 2:
        sum = Decimal(spin_2_coeff_mu * (x**8 + y**8) +
                      spin_2_coeff_nu * (x**6) * (y**2) + spin_2_coeff_nu * (x**2) * (y**6) + 
                      spin_2_coeff_la * (x**4) * (y**4))
        sum/= spin_2_coeff_d

    if spin_S == 3:
        sum = Decimal(spin_3_coeff_mu * (x**8 + y**8) +
                      spin_3_coeff_nu * (x**6) * (y**2) + spin_3_coeff_nu * (x**2) * (y**6) + 
                      spin_3_coeff_la * (x**4) * (y**4))
        sum = Decimal(sum * (x**2 - y**2)**2)
        sum = Decimal(sum / spin_3_coeff_d)

    if spin_S == 4:
        x16 = spin_4_coeff_a * (x**16)
        x14_y2 = -2*spin_4_coeff_b * (x**14 * y**2)



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
    normalized_values = [abs(val / max_value) for val in function_values]
    
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
plt.savefig('plot.png', format='png')