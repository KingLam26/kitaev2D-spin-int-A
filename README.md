# Integer-spin-S Kitaev 2D Honeycomb Model

We consider the arbitary integer-spin-S Kitaev 2D Honeycomb Model [arXiv:1811.05668] in the anisotropic limit using perturbation theory. 

$$
\mathcal{H}_0=-J_z \sum_{\langle i, j\rangle_z} S_i^z S_j^z, \quad V=-J_x \sum_{\langle i, j\rangle_x} S_i^x S_j^x-J_y \sum_{\langle i, j\rangle_y} S_i^y S_j^y .
$$

Explicitly compute the coefficients for each diagram for integer-spin-S kitaev2D case.


## File Structure


The project consists of the following files:

utils.py: Contains utility functions and helper code for processing model parameters, generating combinatorial arrays, and saving results.

core.py: Implements core functionalities, including generating and processing unique permutations with constraints.

main_parallel.py: Orchestrates the parallel processing framework, including producers and consumers for handling permutations across multiple cores.

As of 26/09/24, we have used the code for spin-2, spin-3, spin-4.
Matches Koga's results exactly.
