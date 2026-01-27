import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from scipy.special import erf
try:
    from scipy.integrate import simps
except ImportError:
    from scipy.integrate import simpson as simps

# Procedure for MECP optimization for PCET/PCEnT systems
# We want to find the MECP between the reactant and product diabatic states. 
# The transfering proton is located on the donor for the reactant state and on 
# the acceptor for the product state. The coordinate of other nuclei will be 
# the same for the reactant and the product, i.e., the MECP geometry. 

# To generalize, we want to find the geometry where the energy of the reactant E_R
# and the energy of the product E_P differ by a given amount of energy dE
# This is done by optimizing the objective function with a penalty
# J1 = (E_R + E_P + dE)/2 + coeff * |E_R - E_P - dE|
# or
# J2 = (E_R + E_P + dE)/2 + coeff * (E_R - E_P - dE)^2
# where coeff is a positive number, J depends on the following parameters: 
# the coodinate of the proton on reactant state rp_R, the coordinate of the proton on product state rp_P,
# and the coodinates of other nuclei r_nuc. 

# To take advantage of ASE's optimization module, we create an auxiliary Atoms object, with N+1 atoms
# where N is the total number of atoms in the molecule, the extra atom coresponds to a replica of the
# transfering proton. Thus both rp_R and rp_P are stored as the atomic coordinates in the auxiliary Atoms 
# object. 

# We set the "energy" of the auxiliary Atoms object as the objective function J
# and the "forces" as the negative gradient of E_R, E_P, and J w.r.t rp_R, rp_P, and r_nuc, respectively
# These quantities will be calculated using the MECPPenalty calculator. 

# A pseudo code would look like
# aux_sys = the auxiliary Atoms object
# aux_sys.calc = MECPPenalty(*args, **kwargs)
# apply geometric constraints as needed
# MECPopt = LBFGS(aux_sys, logfile=filename)
# MECPopt.run(fmax=0.05)


class MECPPenalty(Calculator):
    """
    Parameters:
    reactant: an Atoms object for the reactant system (transferring proton on the donor).
    product: an Atoms object for the product system (transferring proton on the acceptor).
    proton_index: (int) index of the transferring proton in the Atoms objects (atomic index starts with 0). 
    reactant_calculator: calculator for the reactant system.
    product_calculator: calculator for the product system.
    penalty_coeff: (float) coefficient for the penalty function, default=1.
    target_dE: (float) targeted energy difference between the reactant and the product at the optimized geometry, dE = E_R - E_P, default = 0.
    opt_state: (string) choose among 'reactant', 'product', and 'average', indicating the energy of which state that the penalty function is added on.
    penalty_function: (string) choose between 'smooth_abs' and 'squared', indicating which type of penalty function is used.

    Optinonal parameters:
    smoothness: (float) control the smoothness of the 'smooth_abs' function when used as the penalty function, smaller value means a smoother function, defult = 1.
    constrain_proton_DA_distance: (Bool) whether to apply constrain on the proton donor-acceptor distance, default = False. 
    donor_index: (int) index of the proton donor (starting with 0). 
    acceptor_index: (int) index of the proton donor (starting with 0). 

    Note that the atomic coodinate of the nuclei other than the transferring proton must be the same 
    in reactant and product. The order of all atoms (including the proton) also must be the same. 

    Currently only implemented for single proton transfer.
    """

    implemented_properties = ['energy', 'forces']

    def __init__(self, reactant, product, proton_index, reactant_calculator, product_calculator, penalty_coeff=1, target_dE=0, state_energies_logfile=None, opt_state='average', penalty_function='squared', constrain_proton_DA_distance=False, **kwargs):
        Calculator.__init__(self, label='MECPPenalty')

        self.reactant = reactant
        self.product = product
        self.proton_index = proton_index
        self.penalty_coeff = penalty_coeff
        self.target_dE = target_dE
        self.reactant_calculator = reactant_calculator
        self.product_calculator = product_calculator
        self.state_energies_logfile = state_energies_logfile
        self.constrain_proton_DA_distance=constrain_proton_DA_distance

        if opt_state not in ['reactant', 'product', 'average']:
            raise ValueError("Choose 'opt_state' in 'reactant', 'product', and 'average'. Default is 'average'. ")
        else:
            self.opt_state = opt_state

        if penalty_function not in ['smooth_abs', 'squared']:
            raise ValueError("Choose 'penalty_function' in 'smooth_abs' and 'squared'. Default is 'squared'. ")
        else:
            self.penalty_function = penalty_function

        if penalty_function == 'smooth_abs':
            if 'smoothness' in kwargs.keys():
                self.smoothness = kwargs['smoothness']
            else:
                self.smoothness = 1

        if constrain_proton_DA_distance:
            if 'donor_index' not in kwargs.keys() or 'acceptor_index' not in kwargs.keys():
                raise NameError("Please provide the index (starting with 0) of the proton donor and acceptor using 'donor_index' and 'acceptor_index'. ")
            else:
                self.donor_index=kwargs['donor_index']
                self.acceptor_index=kwargs['acceptor_index']

        # create the auxiliary Atoms object from reactant and product
        # check atomic coordinates and order
        error_flag = False
        for ii, (atom_R, atom_P) in enumerate(zip(reactant, product)):
            if atom_R.symbol != atom_P.symbol:
                error_flag = True
            if ii != proton_index and (atom_R.position != atom_P.position).all():
                error_flag = True
        if error_flag:
            raise ValueError("The atomic coodinate of the nuclei other than the transferring proton must be the same in reactant and product. The order of all atoms (including the proton) also must be the same.")

        self.aux_system = Atoms()
        self.aux_system += reactant[proton_index]
        self.aux_system += product[proton_index]
        other_nuc = reactant.copy()
        del other_nuc[proton_index]
        self.aux_system += other_nuc


    def calculate(self, atoms, properties=['energy', 'forces'], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        # update the atomic coordinates in self.reactant and self.product using the coordinates in self.atoms
        positions = self.atoms.get_positions()
        current_reactant_positions = np.concatenate((positions[2:self.proton_index+2], positions[[0]], positions[self.proton_index+2:]), axis=0)
        current_product_positions = np.concatenate((positions[2:self.proton_index+2], positions[[1]], positions[self.proton_index+2:]), axis=0)
        
        self.reactant.set_positions(current_reactant_positions)
        self.product.set_positions(current_product_positions)

        self.reactant.calc = self.reactant_calculator
        self.product.calc = self.product_calculator

        self.reactant_calculator.calculate(atoms=self.reactant, properties=['energy', 'forces'])
        self.product_calculator.calculate(atoms=self.product, properties=['energy', 'forces'])
        reactant_energy = self.reactant_calculator.results['energy']
        product_energy = self.product_calculator.results['energy']
        reactant_forces = self.reactant_calculator.results['forces']
        product_forces = self.product_calculator.results['forces']

        if self.state_energies_logfile != None:
            with open(self.state_energies_logfile, 'a') as logfile:
                logfile.write(f'reactant energy = {reactant_energy:16.8f} eV,    product energy = {product_energy:16.8f} eV,    dE = E_reac - E_prod = {reactant_energy-product_energy:8.6f} eV\n')

        if self.penalty_function == 'smooth_abs':
            if self.opt_state == 'reactant':
                # J = E_R + coeff * |E_R - E_P - dE|
                # The term appears in the actural derivetive of J, should be np.sign(reactant_energy - product_energy - self.target_dE)
                # we replace this function by an error function to help convergence
                # This replacement thus changes the actual objective function being optimized to
                # J = E_R + coeff * smooth_abs(E_R - E_P - dE)
                # where smooth_abs(x) = \int_0^x erf(a*x) dx and a is a parameter
                # when a --> \intfy, smooth_abs(x) = |x|
                objective = reactant_energy + self.penalty_coeff * smooth_abs(reactant_energy - product_energy - self.target_dE, a=self.smoothness)
                objective_forces_reac_proton = reactant_forces[[self.proton_index]]
                objective_forces_prod_proton = product_forces[[self.proton_index]]
                objective_forces_other_nuc = np.array([force for i, force in enumerate(reactant_forces \
                                            + self.penalty_coeff*erf(self.smoothness*(reactant_energy - product_energy - self.target_dE)) * (reactant_forces - product_forces)) if i != self.proton_index])
            elif self.opt_state == 'product':
                # J = E_P + coeff * smooth_abs(E_R - E_P - dE)
                # Note that in principle we should define J as J = (E_P + dE) + coeff * smooth_abs(E_R - E_P - dE)
                # The dE term in the first parenthesis is omitted because it's a constant and does not affect the gradient, 
                # and the actural value of J is not physically meaningful
                objective = product_energy + self.penalty_coeff * smooth_abs(reactant_energy - product_energy - self.target_dE, a=self.smoothness)
                objective_forces_reac_proton = reactant_forces[[self.proton_index]]
                objective_forces_prod_proton = product_forces[[self.proton_index]]
                objective_forces_other_nuc = np.array([force for i, force in enumerate(product_forces \
                                            + self.penalty_coeff*erf(self.smoothness*(reactant_energy - product_energy - self.target_dE)) * (reactant_forces - product_forces)) if i != self.proton_index])
            elif self.opt_state == 'average':
                # J = 0.5*(E_R + E_P) + coeff * smooth_abs(E_R - E_P - dE)
                # Note that in principle we should define J as J = 0.5*(E_R + E_P + dE) + coeff * smooth_abs(E_R - E_P - dE)
                # The dE term in the first parenthesis is omitted because it's a constant and does not affect the gradient, 
                # and the actural value of J is not physically meaningful
                objective = 0.5*(reactant_energy + product_energy) + self.penalty_coeff * smooth_abs(reactant_energy - product_energy - self.target_dE, a=self.smoothness)
                objective_forces_reac_proton = reactant_forces[[self.proton_index]]
                objective_forces_prod_proton = product_forces[[self.proton_index]]
                objective_forces_other_nuc = np.array([force for i, force in enumerate(0.5*(reactant_forces + product_forces) \
                                            + self.penalty_coeff*erf(self.smoothness*(reactant_energy - product_energy - self.target_dE)) * (reactant_forces - product_forces)) if i != self.proton_index])
        elif self.penalty_function == 'squared':
            if self.opt_state == 'reactant':
                # J = E_R + coeff * (E_R - E_P - dE)^2
                objective = reactant_energy + self.penalty_coeff * (reactant_energy - product_energy - self.target_dE)**2
                objective_forces_reac_proton = reactant_forces[[self.proton_index]]
                objective_forces_prod_proton = product_forces[[self.proton_index]]
                objective_forces_other_nuc = np.array([force for i, force in enumerate(reactant_forces \
                                            + 2*self.penalty_coeff*(reactant_energy - product_energy - self.target_dE) * (reactant_forces - product_forces)) if i != self.proton_index])
            elif self.opt_state == 'product':
                # J = E_P + coeff * (E_R - E_P - dE)^2
                # Note that in principle we should define J as J = (E_P + dE) + coeff * (E_R - E_P - dE)^2
                # The dE term in the first parenthesis is omitted because it's a constant and does not affect the gradient, 
                # and the actural value of J is not physically meaningful
                objective = product_energy + self.penalty_coeff * (reactant_energy - product_energy - self.target_dE)**2
                objective_forces_reac_proton = reactant_forces[[self.proton_index]]
                objective_forces_prod_proton = product_forces[[self.proton_index]]
                objective_forces_other_nuc = np.array([force for i, force in enumerate(product_forces \
                                            + 2*self.penalty_coeff*(reactant_energy - product_energy - self.target_dE) * (reactant_forces - product_forces)) if i != self.proton_index])
            elif self.opt_state == 'average':
                # J = 0.5*(E_R + E_P) + coeff * (E_R - E_P - dE)^2
                # Note that in principle we should define J as J = 0.5*(E_R + E_P + dE) + coeff * (E_R - E_P - dE)^2
                # The dE term in the first parenthesis is omitted because it's a constant and does not affect the gradient, 
                # and the actural value of J is not physically meaningful
                objective = 0.5*(reactant_energy + product_energy) + self.penalty_coeff * (reactant_energy - product_energy - self.target_dE)**2
                objective_forces_reac_proton = reactant_forces[[self.proton_index]]
                objective_forces_prod_proton = product_forces[[self.proton_index]]
                objective_forces_other_nuc = np.array([force for i, force in enumerate(0.5*(reactant_forces + product_forces) \
                                            + 2*self.penalty_coeff*(reactant_energy - product_energy - self.target_dE) * (reactant_forces - product_forces)) if i != self.proton_index])


        objective_forces = np.concatenate((objective_forces_reac_proton, objective_forces_prod_proton, objective_forces_other_nuc), axis=0)

        if self.constrain_proton_DA_distance:
            rDA = current_reactant_positions[self.acceptor_index] - current_reactant_positions[self.donor_index]
            nDA = rDA/np.linalg.norm(rDA)

            # index of the proton donor and acceptor in the auxiliary system
            idx_D = self.donor_index+2 if self.donor_index < self.proton_index else self.donor_index+1
            idx_A = self.acceptor_index+2 if self.acceptor_index < self.proton_index else self.acceptor_index+1

            objective_forces[idx_D] -= np.inner(objective_forces[idx_D], nDA)*nDA
            objective_forces[idx_A] -= np.inner(objective_forces[idx_A], nDA)*nDA

        self.results['energy'] = objective
        self.results['forces'] = objective_forces


def smooth_abs(x, a=1, ngrid=500):
    if isinstance(x, float):
        xx = np.linspace(0, x, ngrid)
        return simps(erf(a*xx), xx)
    elif isinstance(x, np.ndarray):
        results = np.zeros(len(x))
        for i, xi in enumerate(x):
            xx = np.linspace(0, xi, ngrid)
            results[i] = simps(erf(a*xx), xx)
        return results

