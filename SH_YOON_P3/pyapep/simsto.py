
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
    
    
    def run_mamoen(self, t_max, n_sec = 5):
        t_max_int = np.int32(np.floor(t_max))
        self._n_sec = n_sec
        n_t = t_max_int*n_sec+ 1
        n_comp = self._n_comp
        C_sta = []
        for ii in range(n_comp):
            C_sta.append(self._y_inlet[ii]*self._P_inlet/R_gas/self._T_inlet*1E5)
        t_dom = np.linspace(0,t_max_int, n_t)
        if t_max_int < t_max:
            t_dom = np.concatenate((t_dom, [t_max]))
        
        ## Column design ##
        epsi = self._epsi   # macroscopic void fraction (m^3/m^3)
        A_cr = self._A      # cross-sectional area (m^2)
        L = self._L         # length (m)
        rho_s = self._rho_s # solid density (kg/m^3)
        Cp_gas = self._Cp_g
        Cp_solid = self._Cp_s
        Cv = self._Cv

        ## Adsorbent infor ##
        k_MTC = np.array(self._k_mtc)
        iso_fn = self._iso
        dH_ad = np.array(self._dH)

        ## Feed conditions ##
        P_inlet = self._P_inlet
        T_inlet = self._T_inlet
        y_inlet = np.array(self._y_inlet)

        def massenerbal(y,t):
            C = np.zeros(self._n_comp)
            q = np.zeros(self._n_comp)
            for ii in np.arange(len(q)):
                C[ii] = y[ii]
                q[ii] = y[ii + self._n_comp]
            T = y[-1]         # Temperature (K)
            
            ## Average heat capacity ## 
            Cp_av = np.sum(C*epsi*Cp_gas) + rho_s*(1-epsi)*Cp_solid # J/m^3/K
            Cp_av_inlet = np.sum(Cp_gas*y_inlet)        # J/mol/K
            R_gas = 8.3145  # J/mol/K
            p_gas = C*R_gas*T/1E5 # pressure in (bar)
            P_ov = np.sum(p_gas)

            # Valve operation
            mdot_in = Cv*(P_inlet - P_ov)*P_ov*1E5/R_gas/T # m^3/sec --> mol/sec
            qeq = iso_fn(p_gas, T)
            if len(k_MTC) == 2:
                dqdt = k_MTC[:,0]*(qeq - q) + k_MTC[:,1]*(qeq - q)**2
            else:
                dqdt = k_MTC*(qeq - q)

            dCdt = y_inlet*mdot_in/epsi/A_cr/L - rho_s*(1-epsi)/epsi*dqdt
            dTdt  = Cp_av_inlet*(T_inlet-T)*mdot_in/Cp_av/A_cr/L + rho_s*(1-epsi)/Cp_av*np.sum(dH_ad*dqdt)
            if np.isscalar(dCdt):
                dCdt = [dCdt]
            if np.isscalar(dqdt):
                dqdt = [dqdt]
            dydt = np.concatenate([dCdt, dqdt,[dTdt]])
            
            return dydt
        ## Intiial conditiosn ##
        C_ov = self._P_init/R_gas/self._T_init*1E5
        C_init = np.zeros(self._n_comp)
        for ii in np.arange(self._n_comp):
            C_init[ii] = C_ov*self._y_init[ii]
        y0 = np.concatenate([C_init, self._q_init, [self._T_init]])
        y_outlet = odeint(massenerbal, y0,t_dom)
        self._y = y_outlet
        self._t = t_dom
        return y_outlet, t_dom
    def next_init(self, change_init = True):
        n_comp = self._n_comp
        y = self._y
        C_ov = 0
        q_init = np.zeros(n_comp)
        C_init = np.zeros(n_comp)
        for ii in np.arange(n_comp):
            C_init[ii] = y[-1,ii]
            C_ov = C_ov + C_init[ii]
            q_init[ii] = y[-1,n_comp+ii]
        y_init = C_init/C_ov
        T_init = y[-1,-1]
        P_init = C_ov*R_gas*T_init/1E5
        #q_init = q_init
        if change_init:
            self._P_init = P_init
            self._T_init = T_init
            self._y_init = y_init
            self._q_init = q_init
        return P_init,T_init,y_init,q_init

    def Graph(
        self, index, yaxis_label = None,
        figsize = [7,5],dpi = 85,
        file_name = None, y = None, color = 'k'
        ):
        if y == None:
            y = self._y
        #lstyle = ['-','--','-.',(0,(3,3,1,3,1,3)),':']
        fig,ax = plt.subplots(figsize = figsize, dpi = dpi)
        ax.plot(self._t, y[:,index],
        linewidth = 1.8,
        color = color)
        ax.set_ylabel(yaxis_label,fontsize = 15)
        ax.set_xlabel('time (sec)', fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        plt.grid(ls = ':')
                
        if yaxis_label == None:
            ylab = 'Variable index = {}'.format(index)
            ax.set_ylabel(ylab, fontsize = 15)
        else:
            ax.set_ylabel(yaxis_label, fontsize = 15)
        
        if file_name != None:
            fig.savefig(file_name,bbox_inches='tight')
        return fig, ax

    def Graph_P(
        self, yaxis_label = None,
        figsize = [7,5],dpi = 85, 
        file_name = None, y = None, color = 'k'
        ):
        if y == None:
            y = self._y
        n_comp = self._n_comp
        P = np.zeros_like(self._t)
        for ii in range(n_comp):
            P = P + y[:, ii]*R_gas*y[:, -1]/1E5    # (bar)
        fig,ax = plt.subplots(figsize = figsize, dpi = dpi)
        ax.plot(self._t, P, 
        linewidth = 1.8,
        color = color)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        plt.grid(ls = ':')
        ax.set_ylabel(yaxis_label,fontsize = 15)
        ax.set_xlabel('time (sec)', fontsize = 15)
        #ax.set_title('Pressure', fontsize = 16)
        if yaxis_label == None:
            ylab = 'pressure (bar)'
            ax.set_ylabel(ylab, fontsize = 15)
        else:
            ax.set_ylabel(yaxis_label, fontsize = 15)
        
        if file_name != None:
            fig.savefig(file_name,bbox_inches='tight')

        return fig, ax

    ## Copy the         
    def copy(self):
        import copy
        self_new = copy.deepcopy(self)
        return self_new

def Graph_multi(indx,*argv, yaxis_label = None,
figsize = [7,5],dpi = 85,
file_name = None, color = 'k'):
    y = np.array([argv[0][0][0,indx]])
    t = np.array([argv[0][1][0]])
    t_end = 0 
    for arg in argv:
        y_tmp = arg[0][1:,indx]
        t_tmp = t_end + arg[1][1:]
        y = np.concatenate([y,y_tmp])
        t = np.concatenate([t,t_tmp])
        t_end = t_tmp[-1]
        #y = np.reshape(y,[-1])
        #t = np.reshape(t,[-1])
    fig,ax = plt.subplots(figsize = figsize, dpi = dpi)
    ax.plot(t, y,
    linewidth = 1.8,
    color = color)
    ax.set_ylabel(yaxis_label,fontsize = 15)
    ax.set_xlabel('time (sec)', fontsize = 15)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    plt.grid(ls = ':')
    #fig.legend(fontsize = 14,bbox_to_anchor = legloc)
    
    if yaxis_label == None:
        ylab = 'Variable index = {}'.format(indx)
        ax.set_ylabel(ylab, fontsize = 15)
    else:
        ax.set_ylabel(yaxis_label, fontsize = 15)
    
    if file_name != None:
        fig.savefig(file_name,bbox_inches='tight')

    return fig, ax

def Graph_multi_P(
    *argv, yaxis_label = None,
    figsize = [7,5],dpi = 85,
    file_name = None, color = 'k'):
    n_comp = np.int32((len(argv[0][0][0])-1)/2)
    y = np.array([argv[0][0][0,0:n_comp]])
    T = np.array([argv[0][0][0,-1]])
    t = np.array([argv[0][1][0]])
    
    t_end = 0 
    for arg in argv:
        y_tmp = arg[0][1:,0:n_comp]
        T_tmp = arg[0][1:,-1]
        t_tmp = t_end + arg[1][1:]
        y = np.concatenate([y,y_tmp])
        t = np.concatenate([t,t_tmp])
        T = np.concatenate([T,T_tmp])
        t_end = t_tmp[-1]
        #y = np.reshape(y,[-1])
        #t = np.reshape(t,[-1])
    P = np.sum(y, axis = 1)*R_gas*T/1E5 # (bar)
    fig,ax = plt.subplots(figsize = figsize, dpi = dpi)
    ax.plot(t, P,
    linewidth = 1.8,
    color = color)
    ax.set_ylabel(yaxis_label,fontsize = 15)
    ax.set_xlabel('time (sec)', fontsize = 15)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    plt.grid(ls = ':')
    #fig.legend(fontsize = 14,bbox_to_anchor = legloc)
    
    if yaxis_label == None:
        ylab = 'Pressure (bar)'
        ax.set_ylabel(ylab, fontsize = 15)
    else:
        ax.set_ylabel(yaxis_label, fontsize = 15)
    if file_name != None:
        fig.savefig(file_name,bbox_inches='tight')
    return fig, ax
    
# %% Test the defined class "tank"

# %% Test "adsorbent_info"
if __name__ == '__main__':
    import pyiast
    import pandas as pd
    # Zeolite 13X
    t1 = tank(5,0.031416,2)
    print(t1)

    par_ch4 = [7.26927417, 0.33068804] # Based on mol/kg vs bar
    par_n2 = [0.62864572, 7.26379457, 1.47727665, 0.04093633] # Based on mol/kg vs bar
    dH_list = [16372.5284,11675.65] ## J/mol
    def Lang(p_in,par):
        qtmp = par[0] * par[1]*p_in/(1+par[1]*p_in)
        return qtmp

    def DSLa(p_in,par):
        qtmp1 = par[0]*par[2]*p_in/(1+par[2]*p_in)
        qtmp2 = par[1]*par[3]*p_in/(1+par[3]*p_in)
        qtmp_return = qtmp1 + qtmp2
        return qtmp_return
    P_tmp = np.linspace(0, 49)
    q_ch4 = Lang(P_tmp,par_ch4)
    #print(di_ch4)
    di_ch4 = {'p':P_tmp,
            'q':q_ch4}
    df_ch4 = pd.DataFrame(di_ch4)
    #print(df_ch4)
    iso0 = pyiast.ModelIsotherm(df_ch4,
                                loading_key='q',pressure_key = 'p',
                                model= 'Langmuir', 
                                param_guess = {'M': par_ch4[0],'K':par_ch4[1]})
    q_n2 = DSLa(P_tmp, par_n2)
    di_n2 = {
        'p':P_tmp,
        'q':q_n2}
    df_n2 = pd.DataFrame(di_n2)
    iso1 = pyiast.ModelIsotherm(
        df_n2,
        pressure_key= 'p', loading_key = 'q',
        model = 'DSLangmuir', 
        param_guess = {
            'M1':par_n2[0],
            'M2':par_n2[1],
            'K1':par_n2[2],
            'K2':par_n2[3]
        }
    )
    print(iso0.params)
    print(iso1.params)
    def Arrh(T,T_ref, dH):
        ret = np.exp(np.abs(dH)/8.3145*(1/T - 1/T_ref))
        return ret

    def iso_mix(P,T,):
        P = np.zeros(2)
        P_norm0 = Arrh(T,180,dH_list[0])*P[0]
        P_norm1 = Arrh(T,180,dH_list[1])*P[1]
        P[0] = P_norm0
        P[1] = P_norm1
        if np.sum(P) > 0:
            x_frac = P/np.sum(P)
        else:
            x_frac = np.zeros(2)
        if x_frac[0] < 0.005:
            q = np.zeros(2)
            q[0] = 0
            q[1] = iso0.loading(P_norm0)
            return q
        elif x_frac[1] < 0.005:
            q = np.zeros(2)
            q[0] = iso1.loading(P_norm1)
            q[1] = 0
            return q
        elif x_frac[0] < 0 and x_frac[1] < 0:
            q = np.zeros(2)
            return q
    t1.adsorbent_info(iso_mix,)
    print(t1)

# %% gas property
    molar_mass = [0.016, 0.028]     # kg/mol 
    t1.gas_prop_info(molar_mass)
    print(t1)


# %% mass transfer
    k_mass_trans = [
        [0.00001, 0.00000001],  # [1st: m/sec, 2nd: (m kg)/(sec mol)]
        [0.00001, 0.00000001]]  # [1st: m/sec, 2nd: (m kg)/(sec mol)]
    Surface_A_per_V = 1.5E9     # (1.5 x 10^9) m^2/m^3
    t1.mass_trans_info(k_mass_trans,Surface_A_per_V)
    print(t1)

# %% thermal information
    dH_list = [16372.5284, 11675.65]# J/mol
    T_ref_list = [180,180]          # K
    Cp_g = np.array([40.63,29.22])  # Gas heat capacity: J/mol/K
    Cp_s = 948                      # Solid heat capacity: J/kg/K
    t1.thermal_info(dH_list,Cp_s, Cp_g) 

    print(t1)
# %% feed flow information
    P_inlet = 5         # (bar)
    T_inlet = 298       # (K)
    y_inlet = [0.5, 0.5]    # ([mol/mol, mol/mol])    
    C_valve = 1E-3      # m^3/sec/bar
    t1.feed_flow_info(P_inlet, T_inlet,y_inlet, C_valve)
    print(t1)
# %% Intial conditions
    P_init = 1          # (bar)
    T_init = 298        # (K)
    y_init = [1,0]
    t1.initialC_info(P_init,T_init, y_init)
    print(t1)

# %% Run
    sim_res1 = t1.run_mamoen(200)

# %% Graph
    # Gas components
    fig1,ax1 = t1.Graph(0,'Component 1', color = 'b')
    fig2,ax2 = t1.Graph(1, 'Component 2', color = 'r')
    # Pressure
    #fig_P,ax_P = t1.Graph_P()
    
    # Two component graphs in a single sheet
    fig_tot,ax_tot = plt.subplots()
    ax_tot.plot(ax1.lines[0].get_xdata(), ax1.lines[0].get_ydata(),
    linewidth = 2,
    label = 'Component 1')
    ax_tot.plot(ax2.lines[0].get_xdata(), ax2.lines[0].get_ydata(),
    linewidth = 2,
    label = 'Component 2')
    plt.legend(fontsize = 14)
    ax_tot.set_ylabel(r'concentration (mol/m$^{3}$)',fontsize = 13)
    ax_tot.set_xlabel('time (sec)', fontsize = 13)
    ax_tot.tick_params(axis = 'both', labelsize = 12)
    plt.show()

# %% Second run from the previous time point (t0 = 200) and entire graph
    t1.feed_flow_info(10,300,y_inlet, C_valve)
    t1.next_init()
    sim_res2 = t1.run_mamoen(200)

    fig_multi, ax_multi = Graph_multi(0, sim_res1,sim_res2,
    yaxis_label = 'Gas concentration 1')
    fig_multi_P, ax_multi_P = Graph_multi_P(sim_res1,sim_res2)
    plt.show()