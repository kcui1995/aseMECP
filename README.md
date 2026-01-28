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
We assume that the users are familiar with ASE. Its documentation can be found [here](https://ase-lib.org/about.html).

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

#### Choosing the objective function
The following parameters are used to define the objective function:

8. `opt_state`: (string) choose among 'reactant', 'product', and 'average', indicating the energy of which state that the penalty function is added on.
9. `penalty_function`: (string) choose between 'smooth_abs' and 'squared', indicating which type of penalty function is used.

We have implemented six different objective functions, they all have the form
```math
J = {\rm system\ energy} + c\times{\rm penalty\ function}
```
Since at the crossing point the energies $`E_{\rm R}`$ and $`E_{\rm P}+\Delta`$ should be equal, we can choose to optimize either one or their average. When set `opt_state = reactant`, we set $`{\rm system\ energy}=E_{\rm R}`$ in $`J`$; when set `opt_state = product`, we set $`{\rm system\ energy}=E_{\rm P}`$ in $`J`$; When set `opt_state = average`, we set $`{\rm system\ energy}=(E_{\rm R}+E_{\rm P})/2`$ in $`J`$. The $`\Delta`$ term is omitted because it's a constant and does not affect the gradient, and the actual value of $`J`$ is not physically meaningful. 

The penalty function must be non-negative and reaches zero when the constraint is met. We have implemented two choices. When set `penalty_function = squared`, we set $`{\rm penalty\ function} = ( E_{{\rm R}} - E_{{\rm P}} -\Delta )^2`$ in $`J`$; when set `penalty_function = smooth_abs`, we set $`{\rm penalty\ function} = {\rm smooth\_abs}( E_{{\rm R}} - E_{{\rm P}} -\Delta )`$, where
```math
{\rm smooth\_abs}(x;a) = \int_{0}^{x} {\rm erf}(ax){\rm d}x
```
This is an approximation to the absolute value of $`x`$, but the curve is smoothed near the origin. $`a`$ is the parameter that controlls the smoothness near the origin. When $`a\rightarrow \infty`$, $`{\rm smooth\_abs}(x;\infty)=|x|`$. The parameter `smoothness` sets the value of $`a`$.

10. `smoothness`: (float) control the smoothness of the 'smooth_abs' function when used as the penalty function, smaller value means a smoother function, defult = 1.

> [!TIP]
> In most cases, the default choice of the objective function (`opt_state = average`, `penalty_function = squared`) works well. However, when the system is far from the crossing point, the large force arising from the quadratic penalty function can cause the optimization to explode. Setting `penalty_function = smooth_abs` may solve this issue. 

#### Constrain proton donor-acceptor distance
One can also constrain the distance between the donor and acceptor distance by setting the following parameters:

11. `constrain_proton_DA_distance`: (Bool) whether to apply constrain on the proton donor-acceptor distance, default = False. 
12. `donor_index`: (int) index of the proton donor (starting with 0). 
13. `acceptor_index`: (int) index of the proton donor (starting with 0). 

#### Example

The following code is an example of setting up the `MECPPenalty` calculator. The quantum chemistry calculation will be performed using Q-Chem at TDDFT/6-31+G** level. 
```python
import numpy as np
from ase.io import read, write
from MECP_w_constraints import MECPPenalty
from ase.calculators.qchem import QChem

# read structures of the reactant and the product from xyz files
reac = read('LES_LEPT_geom_move_H.xyz')
prod = read('LEPT_pl_DA_along_Z_aligned.xyz')

# The following lines set up Q-Chem calculaters for the reactant and the product
reac_calc = QChem(label='calc/LES',
                  method='CIS',
                  exchange='CAM-B3LYP',
                  basis='6-31+G**',
                  CIS_N_Roots=5,
                  CIS_triplets='False',
                  CIS_STATE_DERIV=1,
                  RPA='False',
                  nt=16,
                  save=False)

prod_calc = QChem(label='calc/LEPT',
                  method='CIS',
                  exchange='CAM-B3LYP',
                  basis='6-31+G**',
                  CIS_N_Roots=5,
                  CIS_triplets='False',
                  CIS_STATE_DERIV=2,
                  RPA='False',
                  nt=16,
                  save=False)

# The following lines set up the MECPPenalty calculator, using the default objective function and with constraints applied to the
# proton donor-acceptor distance. The calculation results will be logged to the file 'state_energies.log' 
calc = MECPPenalty(reactant=reac,
                    product=prod,
                    proton_index=1,
                    reactant_calculator=reac_calc,
                    product_calculator=prod_calc,
                    state_energies_logfile='state_energies.log',
                    opt_state='average',
                    penalty_function='squared',
                    penalty_coeff=25,
                    constrain_proton_DA_distance=True,
                    donor_index=0,
                    acceptor_index=2)
```

### II. Optimization
With the `MECPPenalty` calculator set up, one can perform the MECP optimization using ASE's optimization module as usual. The following example uses the LBFGS optimizer, and set the convergence criteria to all forces less than 0.03 eV/angstrom.
```python
from ase.optimize import LBFGS

aux_sys = calc.aux_system
aux_sys.calc = calc
opt = LBFGS(aux_sys, logfile='opt.log', trajectory='aux_opt.traj', damping=0.1)
opt.run(fmax=0.03)
```
