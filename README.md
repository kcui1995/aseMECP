# ase-MECP
This code is used for minimum energy crossing point (MECP) optimization specifically for PCET and PCEnT systems. It interfaces with the atomic simulation environment (ASE) to drive first-principles calculations and update the coordinates. 

## Installation 
To use this module, simply download the code and add the parent folder to your `$PYTHONPATH` variable.
For example
```bash
export PYTHONPATH="[YOUR_INSTALLATION_PATH]/ase-MECP-main":$PYTHONPATH
```
You must have the [atomic simulation environment (ASE)](https://ase-lib.org/about.html) installed to use this module. 

The file `qchem.py` is a locally modified version of the original ASE implementation that interfaces Q-Chem with ASE. These modifications enable running CIS calculations in Q-Chem through ASE. When needed, replace the file at `[YOUR_ASE_INSTALLATION_PATH]/ase/calculators/qchem.py` with this version. 

## Background
In the nonadiabatic PCET and PCEnT theories, the transition occurs at the crossing point of the reactant and product diabatic states. Many quantities such as the electronic coupling and the proton potential curves need to be computed at this crossing point to calculate the rate constant. Practically, this crossing point is often approximated as the averaged geometry of the reactant and product structures. However, for certain applications, obtaining the actual crossing point geometry is essential. This code does this job. 

To get the relevant crossing point geometry in the PCET and PCEnT theories, we want to find the coordinates of the nuclei other than the transferring proton such that the energy of the reactant state with the transferring proton optimized near its donor and the energy of the product state with the transferring proton optimized on its acceptor are equal:  
```math
E_{\rm R}(r_{\rm p,eq}^{\rm R},q^*) = E_{\rm P}(r_{\rm p,eq}^{\rm P},q^*)
```
where $`E_{\rm R}`$ and $`E_{\rm P}`$ are the electronic energy of the reactant and product diabatic state, where the electron/excitation energy is localized on the donor and acceptor, respectively, $`q^*`$ denotes the coordinate of all nuceli other than the transferring proton, and $`r_{\rm p,eq}^{\rm R/P}`$ is the equilibrium position of the transferring proton on the reactant/product diabatic state where all other nuclei are fixed at $`q^*`$. We need to optimize $`q^*`$ with this constraint applied to minimize the energies $`E_{\rm R}`$ and $`E_{\rm P}`$. 

To generalize, we want to find the geometry $`q^*`$ where $`E_{\rm R}`$ and $`E_{\rm P}`$ differ by a given amount of energy $\Delta$,
```math
E_{\rm R}(r_{\rm p,eq}^{\rm R},q^*) = E_{\rm P}(r_{\rm p,eq}^{\rm P},q^*) + \Delta
```
This can be used to align arbitrary proton vibrational energy levels associated with the reactant and the product states. 

In practice, $`r_{\rm p,eq}^{\rm R}`$, $`r_{\rm p,eq}^{\rm P}`$, and $`q^*`$ are all unknown and need to be obtained through geometry optimization. To satisfy the condition stated above, we define the objective function where a penalty function is added to the system's energy. This penalty function causes the objective function to increase if the system deviates from the constraint. One possible choice of the objective functions is
```math
J(r_{\rm p1}, r_{\rm p2}, q) = \frac{E_{\rm R}(r_{\rm p1}, q) + E_{\rm P}(r_{\rm p2}, q) + \Delta}{2} + c\left( E_{{\rm R}}(r_{\rm p1}, q) - E_{{\rm P}}(r_{\rm p2}, q) -\Delta \right)^2
```
Here $`r_{\rm p1}`$, $`r_{\rm p2}`$, and $`q`$ are independent variables. Five other choices of the objective functions are also implemented. The restrained optimization is then achieved by minimizing this objective function with respect to $`q`$. The constant $`c`$ controls the strength of the penalty. Meanwhile, we minimize $E_{\rm R}(r_{\rm p1}, q)$ and $E_{\rm P}(r_{\rm p2}, q)$ with respect to $r_{\rm p1}$ and $r_{\rm p2}$ at a given $q$, respectively, to obtain $`r_{\rm p,eq}^{\rm R}`$ and $`r_{\rm p,eq}^{\rm P}`$. The overall optimization process can be summarized as
```math
    \left.\min_{r_{\rm p1}}E_{\rm R}(r_{\rm p1}, q)\right|_{\text{fix }q},\quad \left.\min_{r_{\rm p2}}E_{\rm P}(r_{\rm p2}, q)\right|_{\text{fix }q}, \quad
    \left.\min_{q}\mathcal{J}(r_{\rm p1}, r_{\rm p2}, q)\right|_{\text{fix }r_{\rm p1},r_{\rm p2}}
```
We update the coordinates $`\{r_{\rm p1}, r_{\rm p2}, q\}`$ according to the gradient 
```math
\left\{\frac{\partial E_{\rm R}}{\partial r_{\rm p1}}, \frac{\partial E_{\rm P}}{\partial r_{\rm p2}}, \frac{\partial J}{\partial q}\right\}
```
To take advantage of ASE's optimization module, we create an auxiliary `Atoms` object, with $`N+1`$ atoms where $`N`$ is the total number of atoms in the molecule, the extra atom coresponds to a replica of the transfering proton. Thus both $`r_{\rm p1}`$ and $`r_{\rm p2}`$ are stored as the atomic coordinates in the auxiliary `Atoms` object. 
We set the "energy" of the auxiliary Atoms object as the objective function $`J`$ and the "forces" as $`\left\{\frac{\partial E_{\rm R}}{\partial r_{\rm p1}}, \frac{\partial E_{\rm P}}{\partial r_{\rm p2}}, \frac{\partial J}{\partial q}\right\}`$. These quantities will be calculated using the `MECPPenalty` calculator. Currently it is only implemented for single proton transfer.


## Documentation
We assume that the users are familiar with ASE, its documentation can be found [here](https://ase-lib.org/about.html).

### I. Initialization
 
To perform an MECP optimization, we first set up an `MECPPenalty` object. The required parameters are:

1. `reactant`: an `Atoms` object for the reactant system (transferring proton on the donor).
2. `product`: an `Atoms` object for the product system (transferring proton on the acceptor).
3. `proton_index`: (int) index of the transferring proton in the `Atoms` objects (atomic index starts with 0). 

> [!NOTE]
> The atomic coodinate of the nuclei other than the transferring proton must be the same in `reactant` and `product`. The order of all atoms (including the proton) also must be the same. 

4. `reactant_calculator`: an ASE `calculator` object for the reactant system.
5. `product_calculator`: an ASE `calculator` object for the product system.
6. `penalty_coeff`: (float) coefficient for the penalty function, default = 1.
7. `target_dE`: (float) targeted energy difference between the reactant and the product at the optimized geometry $`\Delta`$, default = 0.
8. `opt_state`: (string) choose among 'reactant', 'product', and 'average', indicating the energy of which state that the penalty function is added on.
9. `penalty_function`: (string) choose between 'smooth_abs' and 'squared', indicating which type of penalty function is used.



Optinonal parameters:

10. `smoothness`: (float) control the smoothness of the 'smooth_abs' function when used as the penalty function, smaller value means a smoother function, defult = 1.
11. `constrain_proton_DA_distance`: (Bool) whether to apply constrain on the proton donor-acceptor distance, default = False. 

When `constrain_proton_DA_distance` is set to `True`, two other parameters are required:

12. `donor_index`: (int) index of the proton donor (starting with 0). 
13. `acceptor_index`: (int) index of the proton donor (starting with 0). 


    
