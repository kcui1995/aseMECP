import numpy as np
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import SCFError
import ase.units


class QChem(FileIOCalculator):
    """
    QChem calculator
    Updated by Kai Cui at Princeton University to support Q-Chem 6 and CIS calculations
    """
    name = 'QChem'

    implemented_properties = ['energy', 'forces']
    command = 'qchem PREFIX.inp PREFIX.out'

    # Following the minimal requirements given in
    default_parameters = {'method': 'hf',
                          'basis': '6-31G*',
                          'jobtype': None,
                          'charge': 0}

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='qchem', scratch=None, np=1, nt=1, pbs=False, save=False,
                 basisfile=None, ecpfile=None, atoms=None, **kwargs):
        """
        The scratch directory, number of processor and threads as well as a few
        other command line options can be set using the arguments explained
        below. The remaining kwargs are copied as options to the input file.
        The calculator will convert these options to upper case
        (Q-Chem standard) when writing the input file.

        scratch: str
            path of the scratch directory
        np: int
            number of processors for the -np command line flag
        nt: int
            number of threads for the -nt command line flag
        pbs: boolean
            command line flag for pbs scheduler (see Q-Chem manual)
        basisfile: str
            path to file containing the basis. Use in combination with
            basis='gen' keyword argument.
        ecpfile: str
            path to file containing the effective core potential. Use in
            combination with ecp='gen' keyword argument.
        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        
        # use new file names for each calculation to avoid overwrite
        self.cnt = 0
        self.label_base = self.label

        # Augment the command by various flags
        if pbs:
            self.command = 'qchem -pbs '
        elif save:
            self.command = 'qchem -save '
        else:
            self.command = 'qchem '
        if np != 1:
            self.command += '-np %d ' % np
        if nt != 1:
            self.command += '-nt %d ' % nt
        self.command += 'PREFIX.inp PREFIX.out'
        if save:
            self.command += ' PREFIX.save'
        if scratch is not None:
            self.command += ' %s' % scratch

        self.basisfile = basisfile
        self.ecpfile = ecpfile

    def read(self, label):
        raise NotImplementedError

    def read_results(self):
        filename = self.label + '.out'

        # if METHOD = CIS then read the excited state energy instead of the ground state energy
        CIS_energy = False
        method_prm = [prm for prm in self.parameters if prm.upper() == 'METHOD'][0]
        if self.parameters[method_prm].upper() == 'CIS':
            CIS_energy = True
            CIS_state = 1           #if 'CIS_STATE_DERIV' is not set, return the energy of the first excited state
            if 'CIS_STATE_DERIV' in [prm.upper() for prm in self.parameters]:
                CIS_state_prm = [prm for prm in self.parameters if prm.upper() == 'CIS_STATE_DERIV'][0]
                CIS_state = self.parameters[CIS_state_prm]

        with open(filename, 'r') as fileobj:
            lineiter = iter(fileobj)
            for line in lineiter:
                if 'SCF failed to converge' in line:
                    raise SCFError()
                elif 'ERROR: alpha_min' in line:
                    # Even though it is not technically a SCFError:
                    raise SCFError()
                elif not CIS_energy and ' Total energy =' in line:
                    convert = ase.units.Hartree
                    self.results['energy'] = float(line.split()[3]) * convert
                elif CIS_energy and f' Total energy for state{CIS_state:4d}' in line: 
                    convert = ase.units.Hartree
                    self.results['energy'] = float(line.split()[5]) * convert
                elif ' Gradient of SCF Energy' in line:
                    # Read gradient as 3 by N array and transpose at the end
                    gradient = [[] for _ in range(3)]
                    # Skip first line containing atom numbering
                    next(lineiter)
                    while True:
                        # Loop over the three Cartesian coordinates
                        for i in range(3):
                            # Cut off the component numbering and remove
                            # trailing characters ('\n' and stuff)
                            line = next(lineiter)[5:].rstrip()
                            # Cut in chunks of 12 symbols and convert into
                            # strings. This is preferred over string.split() as
                            # the fields may overlap for large gradients
                            gradient[i].extend(list(map(
                                float, [line[i:i + 12]
                                        for i in range(0, len(line), 12)])))

                        # After three force components we expect either a
                        # separator line, which we want to skip, or the end of
                        # the gradient matrix which is characterized by the
                        # line ' Max gradient component'.
                        # Maybe change stopping criterion to be independent of
                        # next line. Eg. if not lineiter.next().startswith(' ')
                        if ' Max gradient component' in next(lineiter):
                            # Minus to convert from gradient to force
                            self.results['forces'] = np.array(gradient).T * (-ase.units.Hartree / ase.units.Bohr)
                            break
                elif 'Gradient of the state energy (including CIS Excitation Energy)' in line:
                    # When 'CIS_STATE_DERIV' is set, the output gradient will be the one for the excited state
                    # Read gradient as 3 by N array and transpose at the end
                    gradient = [[] for _ in range(3)]
                    # Skip first line containing atom numbering
                    next(lineiter)
                    while True:
                        # Loop over the three Cartesian coordinates
                        for i in range(3):
                            # Cut off the component numbering and remove
                            # trailing characters ('\n' and stuff)
                            line = next(lineiter)[5:].rstrip()
                            # Cut in chunks of 12 symbols and convert into
                            # strings. This is preferred over string.split() as
                            # the fields may overlap for large gradients
                            gradient[i].extend(list(map(
                                float, [line[i:i + 12]
                                        for i in range(0, len(line), 12)])))

                        # After three force components we expect either a
                        # separator line, which we want to skip, or the end of
                        # the gradient matrix which is characterized by the
                        # line ' Gradient time'.
                        # Maybe change stopping criterion to be independent of
                        # next line. Eg. if not lineiter.next().startswith(' ')
                        if ' Gradient time' in next(lineiter):
                            # Minus to convert from gradient to force
                            self.results['forces'] = np.array(gradient).T * (-ase.units.Hartree / ase.units.Bohr)
                            break

    def write_input(self, atoms, properties=None, system_changes=None):
        self.label = self.label_base + f'{self.cnt:02d}'
        self.cnt += 1

        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        filename = self.label + '.inp'

        with open(filename, 'w') as fileobj:
            fileobj.write('$comment\n   ASE generated input file\n$end\n\n')

            fileobj.write('$rem\n')
            if self.parameters['jobtype'] is None:
                if 'forces' in properties:
                    fileobj.write('   %-25s   %s\n' % ('JOBTYPE', 'FORCE'))
                else:
                    fileobj.write('   %-25s   %s\n' % ('JOBTYPE', 'SP'))

            for prm in self.parameters:
                if prm not in ['charge', 'multiplicity']:
                    if self.parameters[prm] is not None:
                        if isinstance(self.parameters[prm], str):
                            fileobj.write('   %-25s   %s\n' % (
                                prm.upper(), self.parameters[prm].upper()))
                        elif isinstance(self.parameters[prm], int):
                            fileobj.write('   %-25s   %d\n' % (
                                prm.upper(), self.parameters[prm]))
                        elif isinstance(self.parameters[prm], float):
                            fileobj.write('   %-25s   %f\n' % (
                                prm.upper(), self.parameters[prm]))

            # Not even a parameters as this is an absolute necessity
            fileobj.write('   %-25s   %s\n' % ('SYM_IGNORE', 'TRUE'))
            fileobj.write('$end\n\n')

            fileobj.write('$molecule\n')
            # Following the example set by the gaussian calculator
            if ('multiplicity' not in self.parameters):
                tot_magmom = atoms.get_initial_magnetic_moments().sum()
                mult = tot_magmom + 1
            else:
                mult = self.parameters['multiplicity']
            # Default charge of 0 is defined in default_parameters
            fileobj.write('   %d %d\n' % (self.parameters['charge'], mult))
            for a in atoms:
                fileobj.write('   %s  %f  %f  %f\n' % (a.symbol,
                                                       a.x, a.y, a.z))
            fileobj.write('$end\n\n')

            if self.basisfile is not None:
                with open(self.basisfile, 'r') as f_in:
                    basis = f_in.readlines()
                fileobj.write('$basis\n')
                fileobj.writelines(basis)
                fileobj.write('$end\n\n')

            if self.ecpfile is not None:
                with open(self.ecpfile, 'r') as f_in:
                    ecp = f_in.readlines()
                fileobj.write('$ecp\n')
                fileobj.writelines(ecp)
                fileobj.write('$end\n\n')
