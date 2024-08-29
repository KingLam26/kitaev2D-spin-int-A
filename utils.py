import numpy as np
import itertools


##### helper functions ####
def read_model_params(filename):
    """
    Reads and processes model parameters from text file where each line contains a parameter.
    It then returns:
        (1) spin value
        (2) a dictionary of bp values
        (3) a multiplicative factor

    Parameters:
    ----------
    filename : str
        The name of the file containing the model parameters.

    Returns:
    -------
    spin_S : int
        The value of the 'spin_S' parameter as an integer.
        
    bp_dict : dict
        A dictionary containing 'bp_' parameters, where each key is the parameter name
        and each value is a list of integers derived from the corresponding value in the file.

    multi_factor : int
        The value of the 'multi_factor' parameter as an integer.
    
    Example:
    -------
    Given a file 'model_params.txt' with content:
    
    spin_S:         2
    bp_1:           11113333
    bp_2:           00002222

    multi_factor: 513802240
    
    The function will return:
    
    (2, 513802240, {'bp_1': [1, 1, 1, 1, 3, 3, 3, 3], 'bp_2': [0, 0, 0, 0, 2, 2, 2, 2]})
    """

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

def gen_A_combis_dict(combi_sizes):
    """
    Generate a dictionary of valid combinatorial arrays for given sizes.

    For each size in the list of `combi_sizes`, this function creates a base array 
    with an equal number of +1s and -1s (with possibly one extra +1 or -1 to balance 
    the array length). It then generates permutations of the middle part of this 
    base array and filters them based on a rolling sum condition. The result is a 
    dictionary where the key is the size of the combinatorial array and the value 
    is a NumPy array of valid permutations.

    Parameters:
    combi_sizes (list of int): List of sizes for which combinatorial arrays are to be generated.

    Returns:
    dict: A dictionary where each key is a size from `combi_sizes` and the corresponding 
          value is a NumPy array of a physical combinatorial arrays of that size.
    
    Example:
    >>> gen_A_combis_dict([4, 6])
    4 [[-1  1 -1  1]
        [-1  1  1 -1]
        [ 1 -1  1 -1]
        [ 1 -1 -1  1]
        [-1 -1  1  1]
        [ 1  1 -1 -1]]
    6 [[-1  1 -1  1 -1  1]
        [-1  1  1 -1 -1  1]
        [-1 -1  1  1 -1  1]
        [ 1 -1 -1  1  1 -1]
        ...
        [ 1 -1  1 -1  1 -1]]
    """

    A_combis_dict = {}
    
    for size in combi_sizes:
        # Create  base array with n/2 +1s and n/2 -1s
        half_n = size // 2
        base_array = [-1] + [1] * (half_n - 1) + [-1] * (half_n - 1) + [1]

        # Generate combis for the middle part
        middle_part = base_array[1:-1]
        combis = set(itertools.permutations(middle_part))

        # Filter permutations while generating final permutations
        filtered_combis = []
        for combi in combis:
            full_combi = list((-1,) + combi + (1,))
            roll_sum = np.cumsum(full_combi)
            if np.all(roll_sum <= 0):
                filtered_combis.append(full_combi)
                #filtered_combis.append(full_combi[::-1])

        A_combis_dict[size] = np.array(filtered_combis)

    """
    for key, value in A_combis_dict.items():
        print(key, value)
    """
    return A_combis_dict

def save_all_raw(data_list):
    # this will save every bond_permute and output value (Decimal level precision) to a text file
    with open(save_all_raw_filename, 'w') as f:
        for bond_permute, result in data_list:
            f.write(f"{bond_permute} {result}\n")

def save_final_results(bp_label, start_time_str, end_time_str, runtime_minutes, cores, total_sum):
    with open(save_final_results_filename, 'a') as file:
        file.write(f"bp_label: {bp_label}\n")
        file.write(f"start time: {start_time_str}\n")
        file.write(f"end time: {end_time_str}\n")
        file.write(f"run time: {runtime_minutes:.2f} minutes\n")
        file.write(f"cores: {cores}\n")
        file.write(f"total_sum: {total_sum}\n")


##### model parameters ####
# read in from text file and obtain spin_S and bp_dict
spin_S, multi_factor, bp_dict = read_model_params('input-spin-2.txt')

# fixed parameters for all diagrams for integer spin
y_bonds = [0,3]
x_bonds = [1,2]
bond_site_dict = {0: (1,0), 1: (1,5), 2: (3,2), 3: (3,4)}
num_of_sites = 6
y_bond_sites = [site for y_bond in y_bonds for site in bond_site_dict[y_bond]]
gs_energy = spin_S * spin_S * 5
A_combis_dict = gen_A_combis_dict(combi_sizes=[2,4,6,8])

# parameters on saving results
save_final_results_filename = f"output-spin-{spin_S}.txt"
save_all_raw_filename = ""


##### parallel processing parameters #####
num_lockable_lists = 2
queue_maxsize = 10000
num_permutes_print = 5000
compute_switch = True
print_switch_pp = False

if print_switch_pp:
    print(f"num_lockable_lists: {num_lockable_lists}")
    print(f"queue_maxsize: {queue_maxsize}")
    print(f"num_permutes_print: {num_permutes_print}")
    print(f"compute_switch: {str(compute_switch)}")