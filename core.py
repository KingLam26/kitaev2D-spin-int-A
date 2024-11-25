"""
This modeule performs combinatorial and energy computations for spin systems with specific configurations.
It generates permutations of bonds, computes spin configurations, and evaluates energy factors and coefficients
based on quantum mechanical properties of the system. The calculations involve spin ladder operators, site mappings,
and energy propagation across spin sites.

Key Functions:
1. Permutation Generation:
    - `next_unique_permute`: Generates the next lexicographical permutation of an array in-place.
    - `gen_next_permute`: Iteratively yields all unique permutations of a given array.
    - `num_bond_permutes`: Calculates the number of unique permutations for an array with repeated elements.

2. Spin System Processing:
    - `process_permute`: Processes a single bond permutation to calculate the sum of combination coefficients.
    - `compute_combi_coeff`: Calculates the combination coefficient for a specific bond and site configuration.
    - `combi_propagate_one_step`: Handles spin and energy updates for each step of the spin propagation process.
    - `gen_A_sites_combi`: Generates all possible combinations for A-sites based on given replacement rules.
    - `gen_ini_spin_config`: Produces the initial spin configuration based on input sites and actions.

3. Sign Calculations:
    - `gen_OHS_sign`: Computes the overall sign based on the spin and site configuration.

This script can be adapted for quantum mechanical models, particularly those involving spin interactions
and site-specific computations for energy propagation.
"""


##### import #####
from utils import *
import numpy as np
from collections import Counter
import math


##### functions #####
# generate permutes #
def next_unique_permute(array):

    """
    Generate the next lexicographical permutation of the input array in-place.

    This function modifies the input array to produce the next permutation in lexicographical order.
    If the array is already in the highest possible order, the function returns `False` and does not modify the array.
    Otherwise, it returns `True` after rearranging the array into the next permutation.

    Steps:
    1. Find the largest index `k` such that `array[k] < array[k + 1]`. If no such index exists, the array is the last permutation.
    2. Find the largest index `l` greater than `k` such that `array[k] < array[l]`.
    3. Swap `array[k]` with `array[l]`.
    4. Reverse the order of the elements from `array[k + 1]` to the end of the array.

    Args:
        array (list): The list of elements to be permuted.

    Returns:
        bool: `True` if the next permutation is generated, `False` if the array is already in the last permutation order.
    """

    # Step 1: Find the largest index k such that array[k] < array[k + 1]
    k = len(array) - 2
    while k >= 0 and array[k] >= array[k + 1]:
        k -= 1

    if k == -1:
        # This is the last permutation
        return False

    # Step 2: Find the largest index l greater than k such that array[k] < array[l]
    l = len(array) - 1
    while array[k] >= array[l]:
        l -= 1

    # Step 3: Swap the value of array[k] with that of array[l]
    array[k], array[l] = array[l], array[k]

    # Step 4: Reverse the sequence from array[k + 1] up to the end
    array[k + 1:] = array[k + 1:][::-1]  # Use slicing to reverse segment

    return True

def gen_next_permute(array):
    """
    Generate all unique permutations of the input array in lexicographical order.

    This function generates and yields all unique permutations of the input array. The permutations
    are generated in lexicographical order, starting with the smallest permutation (in sorted order).

    Args:
        array (list): The list of elements to permute.

    Yields:
        list: The next permutation of the input list.

    Example:
        >>> for perm in gen_next_permute([1, 2, 3]):
        >>>     print(perm)
        [1, 2, 3]
        [1, 3, 2]
        [2, 1, 3]
        [2, 3, 1]
        [3, 1, 2]
        [3, 2, 1]
    """
    #array.sort()            # Start with the smallest lexicographical permutation
    yield array[:]          # Yield the first permutation
    while next_unique_permute(array):
        yield array[:]

def num_bond_permutes(base):
    n = len(base)
    counts = Counter(base)
    numerator = math.factorial(n)
    denominator = 1
    for count in counts.values():
        denominator *= math.factorial(count)
    num_permutations = numerator // denominator
    return num_permutations

# process a single permute #
def process_permute(bond_permute):

    bond_site_dict_2 = {0: (1,0), 1: (1,5), 2: (3,2), 3: (3,4)}

    Z_sites = np.array([bond_site_dict_2[bond][0] for bond in bond_permute])
    A_sites = np.array([bond_site_dict_2[bond][1] for bond in bond_permute])
    A_sites_combi = gen_A_sites_combi(A_sites)

    # combi computation
    # removed unique_count = len(np.unique(A_sites))
    # no need to consider (2 ** unique_count) multiplicative factor (A sites)
    # and certainly no need for different starting configurations at the Z-links
    Z_fixed_combi_sum = compute_combi_coeff(bond_permute, A_sites_combi, Z_sites, A_sites)

    return bond_permute, Z_fixed_combi_sum

# process a single permute - support #
def gen_A_sites_combi(A_sites):
    # A sites
    A_site_counts = Counter(A_sites)
    A_combis_array = [A_sites]
    for site, count in A_site_counts.items():
        pos_to_replace = np.where(A_sites == site)[0]
        A_combis_array = site_to_ladder(A_combis_array, pos_to_replace, A_combis_dict[count])

    return A_combis_array

def site_to_ladder(combis_array, pos_to_replace, up_down_combis):
    # Ensure that there are positions to replace and permutations to apply
    if len(pos_to_replace) == 0 or len(up_down_combis) == 0:
        return combis_array  # If nothing to replace, return the original array
    replaced_combis = []
    for combi in combis_array:
        for perm in up_down_combis:
            for start_index in range(len(pos_to_replace) - len(perm) + 1):
                # Copy the original array
                new_combi = combi.copy()
                # Place the permutation in the positions of k
                new_combi[pos_to_replace[start_index:start_index + len(perm)]] = perm
                # Append the new array to the list
            replaced_combis.append(new_combi)
    return replaced_combis

def gen_OHS_sign(bond_permute, A_site_combi):
    y_bond_indices = np.where(np.isin(bond_permute, y_bonds))[0]   # no. of y-down (first part)

    # Z-site
    y_down_count_Z = len(y_bond_indices)

    # A-site
    y_bond_up_down_A = A_site_combi[y_bond_indices]
    y_up_count_A = np.sum(y_bond_up_down_A == 1)
    y_down_count_A = np.sum(y_bond_up_down_A == -1)

    # num y_up and y_down
    y_up_count = y_up_count_A
    y_down_count = y_down_count_A + y_down_count_Z
    y_total_count = y_up_count + y_down_count

    sign_1 = 1 if y_total_count % 4 == 0 else -1
    sign_2 = -1 if y_down_count % 2 == 0 else 1
    return sign_1 * sign_2


def gen_ini_spin_config(A_site_combi, A_sites):
    ini_spin_config = [spin_S] * num_of_sites

    # Z sites are all fixed to be +S
    # A sites
    for site in np.unique(A_sites):
        site_first_pos = np.where(A_sites == site)[0][0]
        if A_site_combi[site_first_pos] == 1:
            ini_spin_config[site] = -2
    return tuple(ini_spin_config)


def compute_combi_coeff(bond_permute, A_site_combis, Z_sites, A_sites):

    Z_fixed_combi_sum = 0
    for A_site_combi in A_site_combis:
        ini_spin_config = gen_ini_spin_config(A_site_combi, A_sites)
        spin_config = list(ini_spin_config).copy()
        #print(spin_config)
        cum_energy_factor, cum_spin_factor = 1, 1
        
        for step in range(len(A_site_combi)):

            spin_config, energy_factor, spin_factor = combi_propagate_one_step(ini_spin_config, spin_config, step, Z_sites, A_sites, A_site_combi)
            cum_spin_factor*= spin_factor
            if step != len(A_site_combi)-1:
                cum_energy_factor*=energy_factor
    
        # compute OHS_sign
        OHS_sign = gen_OHS_sign(bond_permute, A_site_combi)

        # compute combi_coeff
        # minus sign for EHS_sign - see v3
        combi_coeff = (OHS_sign * cum_spin_factor) / cum_energy_factor
        print(cum_energy_factor)

        Z_fixed_combi_sum+=combi_coeff
        #print(spin_config)
        #print(bond_permute, A_site_combi, combi_coeff)

    return Z_fixed_combi_sum


def combi_propagate_one_step(ini_spin_config, spin_config, step, Z_sites, A_sites, A_site_combi):

    # extract sites and present spins
    A_site, Z_site = A_sites[step], Z_sites[step]
    A_ini_spin, Z_ini_spin = spin_config[A_site], spin_config[Z_site]
    A_action = A_site_combi[step]

    # compute spin ladder operator factor (without square root)
    # logic block for spin 2
    if spin_S == 2:
        if (A_ini_spin == 1 and A_action == -1) or (A_ini_spin == -1 and A_action == 1) or (A_ini_spin == 0):
            A_spin_factor = 6
        else:
            A_spin_factor = 4

        if Z_ini_spin == 1 or Z_ini_spin == 0:
            Z_spin_factor = 6
        else:
            Z_spin_factor = 4

    # logic block for spin 3
    if spin_S == 3:
        if (A_ini_spin == 1 and A_action == -1) or (A_ini_spin == -1 and A_action == 1) or A_ini_spin == 0:
            A_spin_factor = 12
        elif (A_ini_spin == 2 and A_action == -1) or (A_ini_spin == -2 and A_action == 1) or (A_ini_spin == 1 and A_action == 1) or (A_ini_spin == -1 and A_action == -1):
            A_spin_factor = 10
        else:
            A_spin_factor = 6

        if Z_ini_spin == 1 or Z_ini_spin == 0:
            Z_spin_factor = 12
        elif Z_ini_spin == 2 or Z_ini_spin == -1:
            Z_spin_factor = 10
        else:
            Z_spin_factor = 6

    # logic block for spin 4
    if spin_S == 4:
        if (A_ini_spin == 1 and A_action == -1) or (A_ini_spin == -1 and A_action == 1) or A_ini_spin == 0:
            A_spin_factor = 20
        elif (A_ini_spin == 2 and A_action == -1) or (A_ini_spin == -2 and A_action == 1) or (A_ini_spin == 1 and A_action == 1) or (A_ini_spin == -1 and A_action == -1):
            A_spin_factor = 18
        elif (A_ini_spin == 3 and A_action == -1) or (A_ini_spin == -3 and A_action == 1) or (A_ini_spin == 2 and A_action == 1) or (A_ini_spin == -2 and A_action == -1):
            A_spin_factor = 14
        else:
            A_spin_factor = 8

        if Z_ini_spin == 1 or Z_ini_spin == 0:
            Z_spin_factor = 20
        elif Z_ini_spin == 2 or Z_ini_spin == -1:
            Z_spin_factor = 18
        elif Z_ini_spin == 3 or Z_ini_spin == -2:
            Z_spin_factor = 14
        else:
            Z_spin_factor = 8

    spin_factor = math.sqrt(A_spin_factor * Z_spin_factor)

    # update spin sites
    spin_config[Z_site]+=-1
    spin_config[A_site]+=int(A_action)
    
    # compute energy gap (always positive)
    excited_energy = 0.0
    for site in [0,2,4,5]:
        excited_energy += spin_config[site] * spin_S
    int_z_bond_energy = spin_config[1] * spin_config[3]
    excited_energy+=int_z_bond_energy
    energy_factor = abs(excited_energy - gs_energy)

    # spin_config is now updated spin_config
    # energy_factor and spin_factor are local, i.e., not cumulative
    return spin_config, energy_factor, spin_factor