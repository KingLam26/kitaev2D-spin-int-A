import os
from decimal import Decimal, getcontext
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi
from scipy.optimize import fsolve

getcontext().prec = 30
spin_S = 1.5
consider_all_combis_switch = False


##### Koga results #####
spin_2_coeff_mu = Decimal('18604521')
spin_2_coeff_nu = Decimal('-82758048')
spin_2_coeff_la = Decimal('129273554')
spin_2_coeff_d = Decimal('5138022400')

spin_3_coeff_mu = Decimal('36052814083126422740')
spin_3_coeff_nu = Decimal('-176028114277347622010')
spin_3_coeff_la = Decimal('287126525350219384887')
spin_3_coeff_d = Decimal('152769160756403896320000')

spin_1_5_coeff_mu = Decimal('3214648723397092084')
spin_1_5_coeff_nu = Decimal('1646995686930432837306')
spin_1_5_coeff_la = Decimal('91522768044989658195')
spin_1_5_coeff_denom = Decimal('3076979551468152422400000')

##### Lam results #####
spin_1_5_coeff_a = Decimal('91893.9942230409066477709971778')
spin_1_5_coeff_b = Decimal('-228749.634355255091820103311867')
spin_1_5_coeff_c = Decimal('3185298.38652854643898375623228')
spin_1_5_coeff_d = Decimal('-1352004.10963680688601971400754')
spin_1_5_coeff_e = Decimal('574735.100907972295494161386568')
spin_1_5_coeff_f = Decimal('535358.042592992637547406398269')
spin_1_5_coeff_g = Decimal('575520.108202348095700865382205')

if consider_all_combis_switch == False:
    spin_1_5_coeff_a /= 4
    spin_1_5_coeff_b /= 8
    spin_1_5_coeff_c /= 64
    spin_1_5_coeff_d /= 32
    spin_1_5_coeff_e /= 16
    spin_1_5_coeff_f /= 16
    spin_1_5_coeff_g /= 16

spin_2_coeff_a = Decimal('18604521')
spin_2_coeff_b = Decimal('92064512')
spin_2_coeff_c = Decimal('-41379024')

spin_3_coeff_a = Decimal('36052814083126422740')
spin_3_coeff_b = Decimal('427101825544440584157')
spin_3_coeff_c = Decimal('124066871221800233745')

spin_4_coeff_a = 1.290249294643517e-05
spin_4_coeff_b = 6.0074698994905354e-05 # computed result has minus sign
spin_4_coeff_c = 9.480000447372665e-05
spin_4_coeff_d = 0.0002798033970269736
spin_4_coeff_e = 0.0004415883619644383 # computed result has minus sign
spin_4_coeff_f = 0.0006969445768607373

"""
print(spin_4_coeff_b / spin_4_coeff_a)
print(spin_4_coeff_c / spin_4_coeff_a)
print(spin_4_coeff_d / spin_4_coeff_a)
print(spin_4_coeff_e / spin_4_coeff_a)
print(spin_4_coeff_f / spin_4_coeff_a)

a = spin_4_coeff_a
b = spin_4_coeff_b
c = spin_4_coeff_c
d = spin_4_coeff_d
e = spin_4_coeff_e
f = spin_4_coeff_f

if consider_all_combis_switch:
    a*=4
    b*=8
    c*=8
    d*=16
    e*=16
    f*=16

print(b / a)
print(c / a)
print(d / a)
print(e / a)
print(f / a)

def compute_z(x):
    z = a*(1+x**16) - 2*b*(x**2 + x**14) + (2*c+d)*(x**4 + x **12) + (-2*b-2*e)*(x**6 + x**10) + (2*d+2*a+f)*(x**8)
    return z

def solve_for_x(x0):
    
    # Use fsolve to find the root
    solution = fsolve(compute_z, x0)
    return solution[0]  # Return the solution

x0 = 0.5
x_solution = solve_for_x(x0)

print(x_solution, compute_z(x_solution))

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
    
    if spin_S == 1.5:
        x10_y2 = spin_1_5_coeff_e * (x**10 * y**2)
        x8_y4 = (2*spin_1_5_coeff_b + 2*spin_1_5_coeff_d) * (x**8 + y**4)
        x6_y6 = (spin_1_5_coeff_a + spin_1_5_coeff_c + 2*spin_1_5_coeff_f + 2*spin_1_5_coeff_g) * (x**6 * y**6)
        x4_y8 = (2*spin_1_5_coeff_b + 2*spin_1_5_coeff_d) * (x**4 * y**8)
        x2_y10 = spin_1_5_coeff_e * (x**2 * y**10)

        sum = x10_y2 + x8_y4 + x6_y6 + x4_y8 + x2_y10

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
    
    if spin_S == 4:
        x16 = (spin_4_coeff_a * 4) * (x**16)
        x14_y2 = (-2*spin_4_coeff_b * 8) * (x**14 * y**2)
        x12_y4 = (2*spin_4_coeff_c * 8 + spin_4_coeff_d * 16) * (x**12 * y**4)
        x10_y6 = (-2*spin_4_coeff_b * 8 - 2*spin_4_coeff_e * 16) * (x**10 * y**6)
        x8_y8 = (2*spin_4_coeff_d * 16 + 2*spin_4_coeff_a * 4 + spin_4_coeff_f * 16) * (x**8 * y**8)

        x6_y10 = (-2*spin_4_coeff_b * 8 - 2*spin_4_coeff_e * 16) * (x**6 * y**10)
        x4_y12 = (2*spin_4_coeff_c * 8 + spin_4_coeff_d * 16) * (x**4 * y**12)
        x2_y14 = (-2*spin_4_coeff_b * 8) * (x**2 * y**14)
        y16 = (spin_4_coeff_a * 4) * (y**16)

        sum = x16 + x14_y2 + x12_y4 + x10_y6 + x8_y8 + x6_y10 + x4_y12 + x2_y14 + y16

    return sum

def J_eff_sum_koga(t):

    x, y = trigo(t)

    if spin_S == 1.5:
        sum = Decimal(spin_1_5_coeff_mu * (x**6) * (y**6) +
                      spin_1_5_coeff_nu * (x**2) * (y**2) * ((x**2) - (y**2))**4 +
                      spin_1_5_coeff_la * (x**4) * (y**4) * ((x**2) - (y**2))**2)
        sum/=spin_1_5_coeff_denom
    
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
        x16 = (spin_4_coeff_a) * (x**16)
        x14_y2 = (-2*spin_4_coeff_b) * (x**14 * y**2)
        x12_y4 = (2*spin_4_coeff_c + spin_4_coeff_d) * (x**12 * y**4)
        x10_y6 = (-2*spin_4_coeff_b - 2*spin_4_coeff_e) * (x**10 * y**6)
        x8_y8 = (2*spin_4_coeff_d + 2*spin_4_coeff_a + spin_4_coeff_f) * (x**8 * y**8)

        x6_y10 = (-2*spin_4_coeff_b - 2*spin_4_coeff_e) * (x**6 * y**10)
        x4_y12 = (2*spin_4_coeff_c + spin_4_coeff_d) * (x**4 * y**12)
        x2_y14 = (-2*spin_4_coeff_b) * (x**2 * y**14)
        y16 = (spin_4_coeff_a) * (y**16)

        sum = x16 + x14_y2 + x12_y4 + x10_y6 + x8_y8 + x6_y10 + x4_y12 + x2_y14 + y16

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