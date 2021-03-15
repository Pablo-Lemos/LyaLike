'''
Theory calculation of the simplified lyman alpha parameters
'''

import numpy as np
from typing import Sequence, Union
from cobaya.theory import Theory
from cobaya.likelihood import Likelihood

class Lya(Theory):
    # Options for Pk
    # Default options can be set globally, and updated from requirements as needed

    kmax: float = 0  # Maximum k (1/Mpc units) for Pk, or zero if not needed
    #nonlinear: bool = False  # whether to get non-linear Pk from CAMB/Class
    #z: Union[Sequence, np.ndarray] = []  # redshift sampling
    #extra_args: dict = {}  # extra (non-parameter) arguments passed to ccl.Cosmology()

    kp_kms: float = 0.009 # The wavenumber at which to define the compressed parameters
    z_kms: float = 3.0

    #_default_z_sampling = np.linspace(0, 5, 100)
    z = [0, z_kms]

    def initialize(self):
        self._var_pairs = set()
        self._required_results = {}
    
    def get_requirements(self):
        return {}

    def fit_ln_polynomial(self,xmin,xmax,x,y,deg=2):
        """ Fit a polynomial on the log of the function, within range"""
        x_fit= (x > xmin) & (x < xmax)
        # We could make these less correlated by better choice of parameters
        poly=np.polyfit(np.log(x[x_fit]), np.log(y[x_fit]), deg=deg)
        return np.poly1d(poly)    

    def must_provide(self, **requirements):
        # requirements is dictionary of things requested by likelihoods
        # Note this may be called more than once

        if 'Lya' not in requirements:
            return {}
        options = requirements.get('Lya') or {}
        if 'methods' in options:
            self._required_results.update(options['methods'])

        self.kmax = max(self.kmax, options.get('kmax', self.kmax))
        
        # Dictionary of the things Lya needs from CAMB/CLASS
        needs = {}

        if self.kmax:
            #self.nonlinear = self.nonlinear or options.get('nonlinear', False)
            self._var_pairs.update(
                set((x, y) for x, y in
                    options.get('vars_pairs', [('delta_nonu', 'delta_nonu')])))

            needs['Pk_grid'] = {
                #'vars_pairs': self._var_pairs or [('delta_nonu', 'delta_nonu')],
                #'nonlinear': (True, False) if self.nonlinear else False,
                'vars_pairs': [('delta_nonu', 'delta_nonu')],
                'nonlinear': False,
                'z': self.z,
                'k_max': self.kmax
            }

        needs['Hubble'] = {'z': self.z}
        #needs['comoving_radial_distance'] = {'z': self.z}

        return needs

    #def get_can_support_params(self):
        # return any nuisance parameters that CCL can support
        #return []

    def calculate(self, state, want_derived=True, **params_values_dict):
        if self.kmax:
            for pair in self._var_pairs:
                # Get the matter power spectrum:
                k_Mpc, _, Pk_lin_Mpc = self.provider.get_Pk_grid(var_pair=pair, nonlinear=False)

                #if self.nonlinear:
                    #k, z, Pk_nonlin = self.provider.get_Pk_grid(var_pair=pair, nonlinear=True)

        # use CAMB results to convert Mpc/h to km/s at pivot redshift
        hubble_z=self.provider.get_Hubble(self.z)
        h = hubble_z[0]/100.0
        k_hMpc = k_Mpc/h
        Pk_lin_hMpc = Pk_lin_Mpc[1] *h**3
        H_z = hubble_z[-1]
        dvdX=H_z/(1+self.z[-1])/h

        k_kms=k_hMpc/dvdX
        P_kms=Pk_lin_hMpc*dvdX**3

        # fit polynomial to log power of linear power around kp_kms
        kmin_kms = 0.5*self.kp_kms
        kmax_kms = 2.0*self.kp_kms
        linP_kms_fit=self.fit_ln_polynomial(kmin_kms/self.kp_kms,kmax_kms/self.kp_kms,k_kms/self.kp_kms,P_kms)

        state['Lya'] = {}
        # amplitude of log linear power at kp_kms
        state['Lya']['ln_A_star'] = linP_kms_fit[0]
        # Delta2_star is the dimensionless amplitude at kp_kms
        state['Lya']['Delta2_star'] = np.exp(linP_kms_fit[0])*self.kp_kms**3/(2*np.pi**2)
        # n_star is the effective slope at kp_kms
        state['Lya']['n_star'] = linP_kms_fit[1]
        # alpha_star is the running at kp_kms (twice the curvature)
        state['Lya']['alpha_star'] = 2.0*linP_kms_fit[2]

    def get_Lya(self):
        """
        Get dictionary of Lya computed quantities.
        results['cosmo'] contains the initialized Lya Cosmology object.
        Other entries are computed by methods passed in as the requirements
        :return: dict of results
        """
        return self._current_state['Lya']

class Tester(Likelihood):
    params = {}#'b_hydro': {"prior": {"min": 0, "max": 1}}}

    def get_requirements(self):
        return {'Lya': {"methods": {'test_method': self.test_method},
                        "kmax": 10}}

    def test_method(self, cosmo):
        return None

    def logp(self, true_DL2=0.35,true_neff=-2.3, DL2_err=0.004, neff_err=0.003, r = 0.55, **pars):
        results = self.provider.get_Lya()
        chi2 = ((results['Delta2_star']-true_DL2)**2/DL2_err**2+(results['n_star']-true_neff)**2/neff_err**2
                -2*r*(results['n_star']-true_neff)*(results['Delta2_star']-true_DL2)/DL2_err/neff_err)/(1-r*r)
        return -0.5*chi2
