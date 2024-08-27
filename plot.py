import os
from decimal import Decimal, getcontext
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi

getcontext().prec = 30


##### extract results #####
# Folder path is the current directory where the script is located
folder_path = os.path.dirname(os.path.abspath(__file__))
results_folder_path = os.path.join(folder_path, 'Run 1 - results')

# List of bp indices to skip
skip_indices = {10, 11, 12}

# Dictionary to store the bp index and corresponding decimal value
bp_values = {}

# Loop through bp_1 to bp_19
for i in range(1, 20):
    if i in skip_indices:
        continue
    
    # Construct the filename
    file_name = f'results-2-bp_{i}.txt'
    file_path = os.path.join(results_folder_path, file_name)
    
    # Ensure the file exists before opening
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            # Read all lines in the file
            lines = file.readlines()
            # Extract the 6th line
            if len(lines) >= 6:
                line = lines[5].strip()
                # Extract the number after "Total Sum: "
                if line.startswith("Total Sum:"):
                    value_str = line.split("Total Sum:")[1].strip()
                    # Convert the extracted string to a Decimal
                    decimal_value = Decimal(value_str)
                    # Store in the dictionary with the bp index
                    bp_values[f'bp_{i}'] = decimal_value
            else:
                print(f"File {file_name} does not have enough lines.")
    else:
        print(f"File {file_name} not found.")

# Print or use the extracted values
for key, value in bp_values.items():
    print(f"{key}: {value}")


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

    # Return the negative of the total sum as the coupling constant
    return -total_sum

def J_eff_sum_koga(t):

    x, y = trigo(t)

    d = Decimal('3076979551468152422400000')
    a = Decimal('3214648723397092084')
    b = Decimal('1646995686930432837306')
    c = Decimal('91522768044989658195')

    # Calculate the expression using the given formula
    sum = Decimal(
        a * (x**6) * (y**6) +
        b * (x**2) * (y**2) * ((x**2) - (y**2))**4 +
        c * (x**4) * (y**4) * ((x**2) - (y**2))**2
    )
    
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
plt.plot(t_values / pi, normalized_values_lam, label='lam')
plt.plot(t_values / pi, normalized_values_koga, label='koga')
plt.xlabel('t (radians) / pi')
plt.ylabel('Nomalized J_eff_sum')
plt.title('Nomalized J_eff_sum against theta')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.savefig('normalized_coupling_constant.png', format='png')