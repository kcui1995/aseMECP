# ase-MECP
This code is used for minimum energy crossing point (MECP) optimization specifically for PCET and PCEnT systems. It interfaces with the [atomic simulation environment (ASE)](https://ase-lib.org/about.html) to drive first-principles calculations and update the coordinates. 

## Background
In the nonadiabatic PCET and PCEnT theories, the transition occurs at the crossing point of the reactant and product diabatic states. Many quantities such as the electronic coupling and the proton potential curves need to be computed at this crossing point to calculate the rate constant. Practically, this crossing point is often approximated as the averaged geometry of the reactant and product structures. However, for certain applications such as calculating the indirect coupling in PCEnT reactions, obtaining the actual crossing point geometry is essential. This code does this job. 

In the PCET and PCEnT theories, the transferring proton is treated on the same footing as the electrons. We want to find the coordinates of the nuclei other than the transferring proton such that the energy of the reactant state with the transferring proton optimized near its donor and the energy of the product state with the transferring proton optimized on its acceptor are equal:  
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


## Installation 
To use this module, simply download the code and add the parent folder to your `$PYTHONPATH` variable.
For example
```bash
export PYTHONPATH="[YOUR_INSTALLATION_PATH]/ase-MECP-main":$PYTHONPATH
```
You must have the [atomic simulation environment (ASE)](https://ase-lib.org/about.html) installed to use this module. 

The file `qchem.py` is a locally modified version of the original ASE implementation that interfaces Q-Chem with ASE. These modifications enable running CIS calculations in Q-Chem through ASE. When needed, replace the file at `[YOUR_ASE_INSTALLATION_PATH]/ase/calculators/qchem.py` with this version. 

## Documentation

### I. Initialization
