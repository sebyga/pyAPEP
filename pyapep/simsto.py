
# %% Importing
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# %% Global variable
R_gas = 8.3145      # 8.3145 J/mol/K

# %% Class "tank"

class tank:
    def __init__(self, L, A, n_comp, E_balance = True):
        self._L = L         # length (m)
        self._A = A         # Cross-sectional area (m^2)
        self._n_comp=n_comp # number of gas component
        self._required = {'Design':True,
        'adsorbent_info':False,
        'gas_prop_info': False,
        'mass_trans_info': False,
        }
        if E_balance:
            self._required['thermal_info'] = False
        self._required['feed_flow_info'] = False
        self._required['initialC_info'] = False        
    def __str__(self):
        str_return = '[[Current information included here]] \n'
        for kk in self._required.keys():
            str_return = str_return + '{0:16s}'.format(kk)
            if type(self._required[kk]) == type('  '):
                str_return = str_return+ ': ' + self._required[kk] + '\n'
            elif self._required[kk]:
                str_return = str_return + ': True\n'
            else:
                str_return = str_return + ': False\n'
        return str_return
    def adsorbent_info(self, iso_fn, epsi = 0.3, rho_s = 1000,P_test_range=[0,10], T_test_range = [273,373]):
        T_test = np.linspace(T_test_range[0], T_test_range[1],6)
        p_test = np.zeros([self._n_comp, 6])
        for ii in range(self._n_comp):
            p_tmp = P_test_range[0] + np.random.random(6)*(P_test_range[1] - P_test_range[0])
            p_test[ii,:]=p_tmp
        try:
            for ii,TT in zip(np.arange(6),T_test):
                iso_test = iso_fn(p_test[:,ii], TT)
            if len(iso_test) != self._n_comp:
                print('Output should be a list/narray including {} narray!'.format(self._n_comp))
            else:
                self._iso = iso_fn
                self._rho_s = rho_s
                self._epsi = epsi
                self._required['adsorbent_info'] = True
        except:
            print('You have problem in iso_fn')
            print('Input should be ( [p1_array,p2_array, ...] and T_array )')
            print('Output should be a list/narray including {} narray!'.format(self._n_comp))
    def gas_prop_info(self, Mass_molar):
        stack_true = 0
        if len(Mass_molar) == self._n_comp:
            stack_true = stack_true + 1
        else:
            print('The input variable should be a list/narray with shape ({0:d}, ).'.format(self._n_comp))
        if stack_true == 1:
            self._M_m = Mass_molar
            self._required['gas_prop_info'] = True
    def mass_trans_info(self, k_mass_transfer, a_specific_surf):
        stack_true = 0
        if len(k_mass_transfer) == self._n_comp:
            if np.isscalar(k_mass_transfer[0]):
                order = 1
                self._order_MTC = 1
            else:
                order = 2
                self._order_MTC = 2
            stack_true = stack_true + 1
        else:
            print('The input variable should be a list/narray with shape ({0:d}, ).'.format(self._n_comp))
        if stack_true == 1:
            self._k_mtc = k_mass_transfer
            self._a_surf = a_specific_surf
            self._required['mass_trans_info'] = True

    def feed_flow_info(self, P_inlet, T_inlet, y_inlet, Cv):
        true_stack = 0
        if np.isscalar(P_inlet):
            true_stack = true_stack + 1
        else:
            print("P_inlet should be scalar !!")
        if np.isscalar(T_inlet):
            true_stack = true_stack + 1
        else:
            print("T_inlet should be scalar !!")
        if np.isscalar(Cv):
            true_stack = true_stack + 1
        else:
            print("Cv (valve constant: m^3/sec/bar) should be scalar !!")
        if len(y_inlet) == self._n_comp:
            true_stack = true_stack + 1
        else:
            print("y_in should be [{0:1d},] narray/list !!".format(self._n_comp))            
        if true_stack == 4:
            self._P_inlet = P_inlet
            self._T_inlet = T_inlet
            self._y_inlet = y_inlet
            self._Cv   = Cv
            self._required['feed_flow_info'] = True

    def initialC_info(self, P,T,y,q = None):
        if q == None:
            try:
                q = self._iso(P*np.array(y),T)
            except:
                print('Isotherm model is inappropriate! First use "adsorbent_info."')
                print('Or assign the initial uptake (mol/kg) ')
                return
        if np.isscalar(P) == False:
            print('P should be a scalar.')
            return
        if np.isscalar(T) == False:
            print('P should be a scalar.')
            return
        if len(y) != self._n_comp:
            print('y (gas composition) should be a ({0:1d}) list/narray.'.format(self._n_comp))
            return
        self._P_init = P
        self._T_init = T
        self._y_init = y
        self._q_init = q
        self._required['initialC_info'] = True

    def thermal_info(self, dH_adsorption,
                     Cp_solid, Cp_gas):
        stack_true = 0
        n_comp = self._n_comp
        if len(dH_adsorption) != n_comp:
            print('dH_adsorption should be ({0:d},) list/narray.'.format(n_comp))
        else:
            stack_true = stack_true + 1
        if len(Cp_gas) != n_comp:
            print('Cp_gas should be ({0:d},) list/narray.'.format(n_comp))            
        else:
            stack_true = stack_true + 1
        if stack_true == 2:
            self._dH = dH_adsorption
            self._Cp_s = Cp_solid
            self._Cp_g = Cp_gas
            self._required['thermal_info'] = True
    
    