import os
from decimal import Decimal, getcontext
from math import cos, sin, pi
import matplotlib.pyplot as plt

##### setting parameters #####
getcontext().prec = 30
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

##### helper functions ####
def read_results(filename):
    koga_dict, lam_dict, lam_combi_dict = {}, {}, {}
    with open(filename, 'r') as file:
        for line in file:

            # Strip any leading/trailing whitespace from the line
            line = line.strip()

            # process lines starting with "koga" or "lam"
            if line.startswith('koga') or line.startswith('lam'):
                # Split the line into key and value parts
                key_part, value_str = line.split(':')
                value_str = value_str.strip()

                # Extract the main key (1.5, 2, 3, etc.)
                # e.g., from "koga_1_5_mu", extract "1.5" as main key
                if len(key_part.split('_')) == 4:
                    key_main = key_part.split('_')[1] + '.' + key_part.split('_')[-2]
                else:
                    key_main = key_part.split('_')[1]

                # Extract sub-key (a, b, c, d, e, f, g / mu, nu, la, denom)
                key_sub = key_part.split('_')[-1]

                # Check if there's a comma in the value string
                if ',' in value_str:
                    value_decimal_str, integer_str = value_str.split(',')
                    value_decimal = Decimal(value_decimal_str.strip())
                    integer_value = int(integer_str.strip())
                else:
                    # If there's no integer part after the comma, assume it to be 1
                    value_decimal = Decimal(value_str.strip())

            if line.startswith('koga'):
                # Initialize dictionary for the main key if it doesn't exist
                if key_main not in koga_dict:
                    koga_dict[key_main] = []
                # Assign the decimal value to the appropriate sub-key
                koga_dict[key_main].append(value_decimal)

            elif line.startswith('lam'):
                # Initialize the dictionary for the main key if it doesn't exist
                if key_main not in lam_dict:
                    lam_dict[key_main] = []
                
                if key_main not in lam_combi_dict:
                    lam_combi_dict[key_main] = []                

                # Assign the Decimal and Decimal * value_integer to the appropriate sub-key
                lam_dict[key_main].append(Decimal(value_decimal))
                lam_combi_dict[key_main].append(Decimal(value_decimal*integer_value))

    # convert all lists to tuples in all three dicts
    koga_dict = {k: tuple(v) if isinstance(v, list) else v for k, v in koga_dict.items()}
    lam_dict = {k: tuple(v) if isinstance(v, list) else v for k, v in lam_dict.items()}
    lam_combi_dict = {k: tuple(v) if isinstance(v, list) else v for k, v in lam_combi_dict.items()}

    return koga_dict, lam_dict, lam_combi_dict

def trigo(t):
    x,y = Decimal(cos(t)), Decimal(sin(t))
    return x,y

def J_eff_sum_lam(t, coeff_tuple, spin_S):
    
    x, y = trigo(t)

    # Calculate the sum of terms
    
    if spin_S == '1.5':
        a,b,c,d,e,f,g = coeff_tuple

        x10_y2 = Decimal(e * ((x**10) * (y**2)))
        x8_y4 = Decimal((-2*b - 2*d) * ((x**8) * (y**4)))
        x6_y6 = Decimal((a + c + 2*f + 2*g) * ((x**6) * (y**6)))
        x4_y8 = Decimal((-2*b - 2*d) * ((x**4) * (y**8)))
        x2_y10 = Decimal(e * ((x**2) * (y**10)))

        sum = Decimal(x10_y2 + x8_y4 + x6_y6 + x4_y8 + x2_y10)

    elif spin_S == '2':
        a,b,c = coeff_tuple

        x8 = Decimal(a * (x**8))
        x6_y2 = Decimal(-2*c * ((x**6) * (y**2)))
        x4_y4 = Decimal((2*a + b) * ((x**4) * (y**4)))
        x2_y6 = Decimal(-2*c * ((x**2) * (y**6)))
        y8 = Decimal(a * (y**8))

        sum = Decimal(x8 + y8 + x4_y4 + x6_y2 + x2_y6)

    elif spin_S == '3':
        a,b,c = coeff_tuple

        x12 = Decimal(a * (x**12))
        y12 = Decimal(a * (y**12))
        x6_y6 = Decimal((-2*a - 2*b) * ((x**6) * (y**6)))
        x4_y8 = Decimal((2*c + b) * (x**4 * y**8))
        x8_y4 = Decimal((2*c + b) * (y**4 * x**8))
        x2_y10 = Decimal((-2*c) * (y**10 * x**2))
        x10_y2 = Decimal((-2*c) * (y**2 * x**10))

        sum = Decimal(x12 + y12 + x6_y6 + x4_y8 + x8_y4 + x2_y10 + x10_y2)
    
    elif spin_S == '4':
        a,b,c,d,e,f = coeff_tuple

        x16 = Decimal(a * (x**16))
        x14_y2 = Decimal((-2*b) * (x**14 * y**2))
        x12_y4 = Decimal((2*c + d) * (x**12 * y**4))
        x10_y6 = Decimal((-2*b - 2*e) * (x**10 * y**6))
        x8_y8 = Decimal((2*d + 2*a + f) * (x**8 * y**8))

        x6_y10 = Decimal((-2*b - 2*e) * (x**6 * y**10))
        x4_y12 = Decimal((2*c + d) * (x**4 * y**12))
        x2_y14 = Decimal((-2*b) * (x**2 * y**14))
        y16 = Decimal(a * (y**16))

        sum = Decimal(x16 + x14_y2 + x12_y4 + x10_y6 + x8_y8 + x6_y10 + x4_y12 + x2_y14 + y16)

    return sum

def J_eff_sum_koga(t, coeff_tuple, spin_S):

    x, y = trigo(t)

    if spin_S != '4':
        mu, nu, la, denom = coeff_tuple

    if spin_S == '1.5':
        sum = Decimal(mu * (x**6) * (y**6) +
                      nu * (x**2) * (y**2) * ((x**2) - (y**2))**4 +
                      la * (x**4) * (y**4) * ((x**2) - (y**2))**2)
        sum/=denom
    
    elif spin_S == '2':
        sum = Decimal(mu * (x**8 + y**8) +
                      nu * (x**6) * (y**2) + 
                      nu * (x**2) * (y**6) + 
                      la * (x**4) * (y**4))

    elif spin_S == '3':
        sum = Decimal((mu * (x**8 + y**8) +
                       nu * (x**6) * (y**2) + 
                       nu * (x**2) * (y**6) + 
                       la * (x**4) * (y**4)) * (x**2 - y**2)**2)

    else:
        sum = 1

    return sum

def normalize_J_eff(func, coeff, spin_S, t_values):
    function_values = [Decimal(func(t, coeff, spin_S)) for t in t_values]
    
    # Find maximum magnitude
    max_value = max(abs(val) for val in function_values)
    
    # Normalize the function values
    normalized_J_eff = [abs(val / max_value) for val in function_values]
    
    return normalized_J_eff

def plot_koga_lam(t_values, norm_J_eff_lam, norm_J_eff_lam_combi, norm_J_eff_koga, spin_S):

    plt.figure(figsize=(10, 6))
    plt.plot(t_values / pi, norm_J_eff_lam, label='lam', linewidth = 3)
    plt.plot(t_values / pi, norm_J_eff_lam_combi, label='lam_combi')
    plt.plot(t_values / pi, norm_J_eff_koga, label='koga', linewidth = 1)
    plt.xlabel(r'$t$ (rad) / $pi$')
    plt.ylabel(r'Nomalized $J_{eff}$')
    plt.title(r'Nomalized $J_{eff}$ against $t$' + f': spin-{spin_S}', pad = 10)
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'plot-spin-{spin_S}.png', format='png', dpi=600)
    plt.close()

"""
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