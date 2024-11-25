# Integer-spin-S Kitaev 2D Honeycomb Model

We consider the arbitary integer-spin-S Kitaev 2D Honeycomb Model defined by 

$$H_0=-J_z \sum_{\langle i, j\rangle_z} S_i^z S_j^z, \quad V=-J_x \sum_{\langle i, j\rangle_x} S_i^x S_j^x-J_y \sum_{\langle i, j\rangle_y} S_i^y S_j^y,$$

in the anisotropic limit. The effective Hamiltonian as per perturbation theory is given by

$$H_{\text {eff }}=\sum_{n=1}^{\infty} P V\left(\frac{\mathbb{I}-P}{E_0-\mathcal{H}_0} V\right)^{n-1} P,$$

where $n$ refers to the order in the expansion of $H_{\text {eff }}$, and $P$ is a projection operator onto the subspace spanned by the ground states of $H_0$. In the case of integer spins, the lowest-order ($4S^{\textrm{th}}$) non-trivial effective Hamitonian is given by 

$$H_{\text {eff }} = J_{\text {eff }} \sum_i\sigma_i^x,$$

where the sum is taken over all $z$ bonds in the 2D honeycomb lattice, and $\sigma_i^x$ is a pseudo-half-spin with two degrees-of-freedom. The goal of this repository is to compute the effective coupling constant $J_{\text {eff }}$ for arbitary spin. The set of all allowed perturbative permutation processes contributing non-trivially to $H_{\text {eff }}$ can be organised into a number of diagrams as introduced in [arXiv:1811.05668](
https://doi.org/10.1103/PhysRevB.99.104408). Our code takes a certain diagram and spin-S as input and returns the diagram coefficient given by the sum of all the contributions it contains. As the number of permutation processes to be considered grows by ~$\mathcal{O}(10^3)$ with each spin increment, the problem quickly becomes computationally intractable. Nevertheless, we provide explicit numerical coefficients for up to spin-4, with exact replications of spin-2 and spin-3 as performed in [arXiv:1811.05668](
https://doi.org/10.1103/PhysRevB.99.104408). With further code optimization (e.g., ``Numba`` and vectorization), we expect spin-5 and spin-6 to be tractable with a large cluster.

## File Structure
The project consists of the following files:
1. ``utils.py``: Contains utility functions and helper code for reading in model parameters, generating combinatorial arrays, fine-tuning parallel processing parameters, and saving results.
4. ``core.py``: Implements core functionalities, including generating and processing unique permutations with constraints.
5. ``main_parallel.py``: Orchestrates the parallel processing framework, including producers and consumers for handling permutations across multiple cores.

As of 26/09/24, we have used the code for spin-2, spin-3, spin-4.
Matches Koga's results exactly.
