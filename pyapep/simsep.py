
# %% Import packages
import numpy as np
#from numpy.lib.function_base import _parse_input_dimensions
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
 
# %% Global varaiebls
R_gas = 8.3145      # 8.3145 J/mol/K
 
# %% column (bed) class
def Ergun(C_array,T_array, M_molar, mu_vis, D_particle,epsi_void,d,dd,d_fo, N):
    rho_g = np.zeros(N)
    for c, m in zip(C_array, M_molar):
        rho_g = rho_g + m*c
    P = np.zeros(N)
    for cc in C_array: 
        P = P + cc*R_gas*T_array
    dP = d@P
    
    Vs_Vg = (1-epsi_void)/epsi_void
    A =1.75*rho_g/D_particle*Vs_Vg
    B = 150*mu_vis/D_particle**2*Vs_Vg**2
    C = dP
 
    ind_posi = B**2-4*A*C >= 0
    ind_nega = ind_posi == False
    v_pos = (-B[ind_posi]+ np.sqrt(B[ind_posi]**2-4*A[ind_posi]*C[ind_posi]))/(2*A[ind_posi])
    v_neg = (B[ind_nega] - np.sqrt(B[ind_nega]**2+4*A[ind_nega]*C[ind_nega]))/(2*A[ind_nega])
    
    v_return = np.zeros(N)
    v_return[ind_posi] = v_pos
    v_return[ind_nega] = v_neg
    
    dv_return = d_fo@v_return
    return v_return, dv_return

def Ergun_test(dP,M_molar, mu_vis, D_particle,epsi_void):
    rho_g = 40*M_molar
    Vs_Vg = (1-epsi_void)/epsi_void
    A =1.75*rho_g/D_particle*Vs_Vg*np.ones_like(dP)
    B = 150*mu_vis/D_particle**2*Vs_Vg**2*np.ones_like(dP)
    C = dP
 
    ind_posi = B**2-4*A*C >= 0
    ind_nega = ind_posi == False
    v_pos = (-B[ind_posi]+ np.sqrt(B[ind_posi]**2-4*A[ind_posi]*C[ind_posi]))/(2*A[ind_posi])
    v_neg = (B[ind_nega] - np.sqrt(B[ind_nega]**2+4*A[ind_nega]*C[ind_nega]))/(2*A[ind_nega])
    
    v_return = np.zeros_like(dP)
    v_return[ind_posi] = v_pos
    v_return[ind_nega] = v_neg
    
    return v_return


def change_node_fn(z_raw, y_raw, N_new):
    if isinstance(y_raw,list):
        fn_list = []
        y_return = []
        z_new = np.linspace(z_raw[0], z_raw[-1],N_new)
        for yr in y_raw:
            fn_tmp = interp1d(z_raw,yr)
            y_new_tmp = fn_tmp(z_new)
            y_return.append(y_new_tmp)
    elif len(y_raw.shape) == 1:
        yy = y_raw
        fn_tmp = interp1d(z_raw, yy, kind = 'cubic')
        z_new = np.linspace(z_raw[0], z_raw[-1],N_new)
        y_return = fn_tmp(z_new)
    elif len(y_raw.shape) == 2:
        yy = y_raw[-1,:]
        fn_tmp = interp1d(z_raw, yy, kind = 'cubic')
        z_new = np.linspace(z_raw[0], z_raw[-1],N_new)
        y_return = fn_tmp(z_new)
    else:
        print('Input should be 1d or 2d array.')
        return None
    return y_return

# %% Column class
class column:
    def __init__(self, L, A_cross, n_component, 
                 N_node = 11, E_balance = True):
        self._L = L
        self._A = A_cross
        self._n_comp = n_component
        self._N = N_node
        self._z = np.linspace(0,L,N_node)
        self._required = {'Design':True,
        'adsorbent_info':False,
        'gas_prop_info': False,
        'mass_trans_info': False,}
        if E_balance:
            self._required['thermal_info'] = False
        self._required['boundaryC_info'] = False
        self._required['initialC_info'] = False
        h = L/(N_node-1)
        h_arr = h*np.ones(N_node)
        self._h = h
 
        # FDM backward, 1st deriv
        d0 = np.diag(1/h_arr, k = 0)
        d1 = np.diag(-1/h_arr[1:], k = -1)
        d = d0 + d1
        d[0,:] = 0
        self._d = d
        
        # FDM foward, 1st deriv
        d0_fo = np.diag(-1/h_arr, k = 0)
        d1_fo = np.diag(1/h_arr[1:], k = 1)
        d_fo = d0_fo + d1_fo
        self._d_fo  = d_fo
 
        # FDM centered, 2nd deriv
        dd0 = np.diag(1/h_arr[1:]**2, k = -1)
        dd1 = np.diag(-2/h_arr**2, k = 0)
        dd2 = np.diag(1/h_arr[1:]**2, k = 1)
        dd = dd0 + dd1 + dd2
        dd[0,:]  = 0
        dd[-1,:] = 0
        self._dd = dd
 
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
 
 ### Before running the simulations ###

    def adsorbent_info(self, iso_fn, epsi = 0.3, D_particle = 0.01, rho_s = 1000,P_test_range=[0,10], T_test_range = [273,373]):
        T_test = np.linspace(T_test_range[0], T_test_range[1],self._N)
        p_test = []
        for ii in range(self._n_comp):
            p_tmp = P_test_range[0] + np.random.random(self._N)*(P_test_range[1] - P_test_range[0])
            p_test.append(p_tmp)        
        try:      
            iso_test = iso_fn(p_test, T_test)
            if len(iso_test) != self._n_comp:
                print('Output should be a list/narray including {} narray!'.format(self._n_comp))
            else:
                self._iso = iso_fn
                self._rho_s = rho_s
                self._epsi = epsi
                self._D_p =D_particle
                self._required['adsorbent_info'] = True
        except:
            print('You have problem in iso_fn')
            print('Input should be ( [p1_array,p2_array, ...] and T_array )')
            print('Output should be a list/narray including {} narray!'.format(self._n_comp))
        
        
    def gas_prop_info(self, Mass_molar, mu_viscosity):
        stack_true = 0
        if len(Mass_molar) == self._n_comp:
            stack_true = stack_true + 1
        else:
            print('The input variable should be a list/narray with shape ({0:d}, ).'.format(self._n_comp))
        if len(mu_viscosity) == self._n_comp:
            stack_true = stack_true + 1
        else:
            print('The input variable should be a list/narray with shape ({0:d}, ).'.format(self._n_comp))
        if stack_true == 2:
            self._M_m = Mass_molar
            self._mu = mu_viscosity
            self._required['gas_prop_info'] = True
    def mass_trans_info(self, k_mass_transfer, a_specific_surf, D_dispersion = 1E-8):
        stack_true = 0
        if np.isscalar(D_dispersion):
            D_dispersion = D_dispersion*np.ones(self._n_comp)
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
        if len(D_dispersion) == self._n_comp:
            stack_true = stack_true + 1
        else:
            print('The input variable should be a list/narray with shape ({0:d}, ).'.format(self._n_comp))
        if stack_true == 2:
            self._k_mtc = k_mass_transfer
            self._a_surf = a_specific_surf
            self._D_disp = D_dispersion
            self._required['mass_trans_info'] = True

    
    def thermal_info(self, dH_adsorption,
                     Cp_solid, Cp_gas, h_heat_transfer,
                     k_conduct = 0.0001, h_heat_ambient = 0, T_ambient = 298):
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
        if np.isscalar(h_heat_transfer):
            stack_true = stack_true + 1
        else:
            print('h_heat_transfer should be scalar.')
        if stack_true == 3:
            self._dH = dH_adsorption
            self._Cp_s = Cp_solid
            self._Cp_g = Cp_gas
            self._h_heat = h_heat_transfer
            self._k_cond = k_conduct
            self._h_ambi = h_heat_ambient
            self._T_ambi = T_ambient
            self._required['thermal_info'] = True
 
    def boundaryC_info(self, P_outlet,
    P_inlet,T_inlet,y_inlet, Cv_in=1E-1,Cv_out=1E-3,
      Q_inlet=None, assigned_v_option = True, 
      foward_flow_direction =  True):
        self._Q_varying = False
        self._required['Flow direction'] = 'Foward'
        if foward_flow_direction == False:
            A_flip = np.zeros([self._N,self._N])
            for ii in range(self._N):
                A_flip[ii, -1-ii] = 1
            self._required['Flow direction'] = 'Backward'
            self._A_flip = A_flip
        if Q_inlet == None:
            assigned_v_option = False
        elif np.isscalar(Q_inlet) == False:
            assigned_v_option = True
            t = Q_inlet[0]
            Q = Q_inlet[1]
            f_Q_in = interp1d(t,Q)
            self._fn_Q = f_Q_in
            self._Q_varying = True
        try:
            if len(y_inlet) == self._n_comp:
                self._P_out = P_outlet
                self._P_in = P_inlet
                self._T_in = T_inlet
                self._y_in = y_inlet
                self._Q_in = Q_inlet
                self._Cv_in = Cv_in
                self._Cv_out = Cv_out
                self._const_v = assigned_v_option
                self._required['boundaryC_info'] = True
                if assigned_v_option:
                    self._required['Assigned velocity option'] = True
                else:
                    self._required['Assigned velocity option'] = False  
            else:
                print('The inlet composition should be a list/narray with shape (2, ).')
        except:
            print('The inlet composition should be a list/narray with shape (2, ).')    
 
    def initialC_info(self,P_initial, Tg_initial,Ts_initial, y_initial,q_initial):
        stack_true = 0
        if len(P_initial) != self._N:
            print('P_initial should be of shape ({},)'.format(self._N))
        else:
            stack_true = stack_true + 1
        if len(y_initial) != self._n_comp or len(y_initial[0]) != self._N:
            print('y_initial should be a list including {0} ({1},) array'.format(self._n_comp, self._N))
        else:
            stack_true = stack_true + 1
        if len(q_initial) != self._n_comp or len(q_initial[0]) != self._N:
            print('q_initial should be a list/array including {0} ({1},) array'.format(self._n_comp, self._N))
        else:
            stack_true = stack_true + 1
        if stack_true == 3:
            self._P_init = P_initial
            self._Tg_init = Tg_initial
            self._Ts_init = Ts_initial
            self._y_init = y_initial
            self._q_init = q_initial
            self._required['initialC_info'] = True               

#########################
##### RUN FUNCTIONS #####
#########################

## Run mass & momentum balance equations
    def run_mamo(self, t_max, n_sec = 5, CPUtime_print = False):
        tic = time.time()/60
        t_max_int = np.int32(np.floor(t_max), )
        self._n_sec = n_sec
        n_t = t_max_int*n_sec+ 1
        n_comp = self._n_comp
        C_sta = []
        for ii in range(n_comp):
            C_sta.append(self._y_in[ii]*self._P_in/R_gas/self._T_in*1E5)
        t_dom = np.linspace(0,t_max_int, n_t)
        if t_max_int < t_max:
            t_dom = np.concatenate((t_dom, [t_max]))
        #print(C1_sta)

        # other parmeters
        epsi = self._epsi
        rho_s = self._rho_s

        N = self._N
        def massmomebal(y,t):
            C = []
            q = []
            for ii in range(n_comp):
                C.append(y[ii*N:(ii+1)*N])
                q.append(y[n_comp*N + ii*N : n_comp*N + (ii+1)*N])
 
            # Derivatives
            dC = []
            ddC = []
            C_ov = np.zeros(N)
            P_ov = np.zeros(N)
            P_part = []
            Mu = np.zeros(N)
            T = self._Tg_init
            for ii in range(n_comp):
                dC.append(self._d@C[ii])
                ddC.append(self._dd@C[ii])
                P_part.append(C[ii]*R_gas*T/1E5) # in bar
                C_ov = C_ov + C[ii]
                P_ov = P_ov + C[ii]*R_gas*T
                Mu = C[ii]*self._mu[ii]
            Mu = Mu/C_ov

            # Ergun equation
            v,dv = Ergun(C,T,self._M_m,Mu,self._D_p,epsi,
                         self._d,self._dd,self._d_fo, self._N)
            
            # Solid phase
            qsta = self._iso(P_part, T) # partial pressure in bar
            dqdt = []
            if self._order_MTC == 1:
                for ii in range(n_comp):
                    dqdt_tmp = self._k_mtc[ii]*(qsta[ii] - q[ii])*self._a_surf
                    dqdt.append(dqdt_tmp)
            elif self._order_MTC == 2:
                for ii in range(n_comp):
                    dqdt_tmp = self._k_mtc[ii][0]*(qsta[ii] - q[ii])*self._a_surf + self._k_mtc[ii][1]*(qsta[ii] - q[ii])**2*self._a_surf
                    dqdt.append(dqdt_tmp)
            # Valve equations (v_in and v_out)
            #P_in = (C1_sta + C2_sta)*R_gas*T_gas
            if self._const_v:
                v_in = self._Q_in/epsi/self._A
            else:
                v_in = max(self._Cv_in*(self._P_in - P_ov[0]/1E5), 0 )  # pressure in bar           
            v_out = max(self._Cv_out*(P_ov[-1]/1E5 - self._P_out), 0 )  # pressure in bar
            
            D_dis = self._D_disp
            h = self._h
            # Gas phase
            dCdt = []
            for ii in range(n_comp):
                dCdt_tmp = -v*dC[ii] -C[ii]*dv + D_dis[ii]*ddC[ii] - (1-epsi)/epsi*rho_s*dqdt[ii]
                #dCdt_tmp[0] = +(v_in*C_sta[ii] - v[1]*C[ii][0])/h - (1-epsi)/epsi*rho_s*dqdt[ii][0]
                #dCdt_tmp[-1]= +(v[-1]*C[ii][-2]- v_out*C[ii][-1])/h - (1-epsi)/epsi*rho_s*dqdt[ii][-1]
                dCdt_tmp[0] = v_in*(C_sta[ii]-C[ii][0])/h- (1-epsi)/epsi*rho_s*dqdt[ii][0]
                dCdt_tmp[-1]= +(v_out+v[-1])/2*(C[ii][-2]- C[ii][-1])/h - (1-epsi)/epsi*rho_s*dqdt[ii][-1]
                dCdt.append(dCdt_tmp)
 
            dydt_tmp = dCdt+dqdt
            dydt = np.concatenate(dydt_tmp)
            return dydt
        
        C_init = []
        q_init = []
        for ii in range(n_comp):
            C_tmp = self._y_init[ii]*self._P_init*1E5/R_gas/self._Tg_init
            C_init.append(C_tmp)
            q_tmp = self._q_init[ii]
            q_init.append(q_tmp)
 
        tic = time.time()/60
        if self._required['Flow direction'] == 'Backward':
            for ii in range(n_comp):
                C_init[ii] = self._A_flip@C_init[ii]
                q_init[ii] = self._A_flip@q_init[ii]
        
        y0_tmp = C_init + q_init
        y0 = np.concatenate(y0_tmp)
        
        #RUN
        y_result = odeint(massmomebal,y0,t_dom,)
        
        if self._required['Flow direction'] == 'Backward':
            y_tmp = []
            for ii in range(n_comp*2):
                mat_tmp = y_result[:, ii*N : (ii+1)*N]
                y_tmp.append(mat_tmp@self._A_flip)
            y_flip = np.concatenate(y_tmp, axis = 1)
            y_result = y_flip
        self._y = y_result
        self._t = t_dom
        toc = time.time()/60 - tic
        self._CPU_min = toc
        self._Tg_res = np.ones([len(self._t), 1])@np.reshape(self._Tg_init,[1,-1])

        if CPUtime_print:
            print('Simulation of this step is completed.')
            print('This took {0:9.3f} mins to run. \n'.format(toc))       
        return y_result, self._z, t_dom

    def Q_valve(self, draw_graph = False, y = None):
        N = self._N
        if self._required['Flow direction'] == 'Backward':
            Cv_0 = self._Cv_out
            Cv_L = self._Cv_in
        else:
            Cv_0 = self._Cv_in
            Cv_L = self._Cv_out
            
        if y == None:
            y = self._y
        P_0 = (y[:,0] + y[:,N])*R_gas*self._Tg_res[:,0]/1E5
        P_L = (y[:,N-1] + y[:,2*N-1])*R_gas*self._Tg_res[:,-1]/1E5
        
        if self._required['Flow direction'] == 'Backward':
            if self._required['Assigned velocity option']:
                v_L = self._Q_in*np.ones(len(self._t))/self._A/self._epsi
            else:
                v_L = Cv_L*(self._P_in - P_L)
            v_0 = Cv_0*(P_0 - self._P_out)
        else:
            if self._required['Assigned velocity option']:
                v_0 = self._Q_in*np.ones(len(self._t))/self._A/self._epsi
            else:
                v_0 = Cv_0*(self._P_in - P_0)
            v_L = Cv_L*(P_L - self._P_out)
        Q_0 = v_0*self._A * self._epsi
        Q_0[Q_0 < 0] = 0
        Q_L = v_L*self._A * self._epsi
        Q_L[Q_L < 0] = 0
 
        if draw_graph:
            plt.figure(figsize = [6.5,5], dpi = 90)
            plt.plot(self._t, Q_0,
            label = 'Q at z = 0', linewidth = 2)
            plt.plot(self._t, Q_L, 
            label = 'Q at z = L', linewidth = 2)
            plt.legend(fontsize = 14)
            plt.xlabel('time (sec)', fontsize = 15)
            plt.ylabel('volumetric flowrate (m$^{3}$/sec)', fontsize =15)
            plt.xticks(fontsize = 14)
            plt.yticks(fontsize = 14)
        return Q_0, Q_L
    
    def breakthrough(self, draw_graph = False):
        N = self._N
        n_comp = self._n_comp
        C_end =[]
        C_ov = np.zeros(len(self._t))
        for ii in range(n_comp):
            Cend_tmp = self._y[:,(ii+1)*N-1]
            C_end.append(Cend_tmp)
            C_ov = C_ov + Cend_tmp
        y = []
        fn_y = [] 
        for ii in range(n_comp):
            y_tmp = C_end[ii]/C_ov
            y.append(y_tmp)
            fn_tmp = interp1d(self._t, y_tmp, kind = 'cubic' ) 
            fn_y.append(fn_tmp)
        if draw_graph:
            t_dom = np.linspace(self._t[0], self._t[-1],1001)
            y_again = []
            plt.figure(figsize = [8,5], dpi = 90)
            for ii in range(n_comp):
                y_ag_tmp = fn_y[ii](t_dom)
                y_again.append(y_ag_tmp)
                plt.plot(t_dom, y_ag_tmp,
                         label = 'Component{0:2d}'.format(ii+1),
                         linewidth = 2 )            
            plt.legend(fontsize = 14)
            plt.xlabel('time (sec)', fontsize = 15)
            plt.ylabel('mole fraction (mol/mol)', fontsize =15)
            plt.xticks(fontsize = 14)
            plt.yticks(fontsize = 14)
    
        return fn_y
 
    def Graph(self, every_n_sec, index, 
              loc = [1,1], yaxis_label = None, 
              file_name = None, 
              figsize = [7,5], dpi = 85, y = None,):
        N = self._N
        one_sec = self._n_sec
        n_show = one_sec*every_n_sec
        if y == None:
            y = self._y
        lstyle = ['-','--','-.',(0,(3,3,1,3,1,3)),':']
        fig, ax = plt.subplots(figsize = figsize, dpi = 90)
        cc= 0
        lcolor = 'k'
        for j in range(0,len(self._y), n_show):
            if j <= 1:
                lcolor = 'r'
            elif j >= len(self._y)-n_show:
                lcolor = 'b'
            else:
                lcolor = 'k'
            ax.plot(self._z,self._y[j, index*N:(index+1)*N],
            color = lcolor, linestyle = lstyle[cc%len(lstyle)],
            label = 't = {}'.format(self._t[j]))
            cc = cc + 1
        fig.legend(fontsize = 14,bbox_to_anchor = loc)
        ax.set_xlabel('z-domain (m)', fontsize = 15)
        if yaxis_label == None:
            ylab = 'Variable index = {}'.format(index)
            ax.set_ylabel(ylab, fontsize = 15)
        else:
            ax.set_ylabel(yaxis_label, fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        plt.grid(linestyle = ':')
        if file_name != None:
            fig.savefig(file_name, bbox_inches='tight')
        
        return fig, ax
        
    def Graph_P(self, every_n_sec, loc = [1,1], 
                yaxis_label = 'Pressure (bar)',
                file_name = None, 
                figsize = [7,5], dpi = 85, y = None,):
        N = self._N
        one_sec = self._n_sec
        n_show = one_sec*every_n_sec
        if y == None:
            y = self._y
        lstyle = ['-','--','-.',(0,(3,3,1,3,1,3)),':']
        P = np.zeros(N)
        for ii in range(self._n_comp):
            P = P + self._y[:,(ii)*N:(ii+1)*N]*R_gas*self._Tg_res/1E5
        fig, ax = plt.subplots(figsize = figsize, dpi = dpi)
        cc= 0
        for j in range(0,len(self._y), n_show):
            if j <= 1:
                lcolor = 'r'
            elif j >= len(self._y)-n_show:
                lcolor = 'b'
            else:
                lcolor = 'k'
            ax.plot(self._z,P[j, :],
            color = lcolor, linestyle = lstyle[cc%len(lstyle)],
            label = 't = {}'.format(self._t[j]))
            cc = cc + 1
        fig.legend(fontsize = 14, bbox_to_anchor = loc)
        ax.set_xlabel('z-domain (m)', fontsize = 15)
        ax.set_ylabel(yaxis_label, fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        plt.grid(linestyle = ':')
        if file_name != None:
            fig.savefig(file_name, bbox_inches='tight')
        return fig, ax    


# %% When only this code is run (name == main)
if __name__ == '__main__':
    N = 21
    A_cros = 0.031416
    L = 1
    c1 = column(L,A_cros, n_component = 2,N_node = N)

    ## Adsorbent
    isopar1 = [3.0, 1]
    isopar2 = [1.0, 0.5]
    def iso_fn_test(P,T):
        b1 = isopar1[1]*np.exp(30E3/8.3145*(1/T-1/300))
        b2 = isopar2[1]*np.exp(20E3/8.3145*(1/T-1/300))
        denom = 1 + b1*P[0] + b2*P[1]
        numor0 = isopar1[0]*b1*P[0]
        numor1 = isopar2[0]*b2*P[1]
        q_return = [numor0/denom, numor1/denom]
        return q_return
 
    epsi_test = 0.4         # macroscopic void fraction (m^3/m^3)
    D_particle_dia = 0.01   # particle diameter (m)
    rho_s_test = 1100       # solid density (kg/mol)
    c1.adsorbent_info(iso_fn_test,epsi_test,D_particle_dia, rho_s_test)
 
    M_m_test  = [0.044, 0.028]      ## molar mass    (kg/mol)
    mu_test = [1.47E-5, 1.74E-5]    ## gas viscosity (Pa sec)
    c1.gas_prop_info(M_m_test, mu_test)
 
    ## Mass transfer coefficients
    D_dis_test = [1E-6, 1E-6]   # m^2/sec
    k_MTC = [1.0, 1.0]        # m/sec
    a_surf = 400                # m^2/m^3
    c1.mass_trans_info(k_MTC, a_surf, D_dis_test)

    ## Thermal properties
    Del_H = [30E3, 20E3]    # J/mol
    Cp_s = 935              # J/kg/K
    Cp_g = [37.22, 29.15]   # J/mol/K
    h_heat = 100            # J/sec/m^2/K
    c1.thermal_info(Del_H,Cp_s,Cp_g,h_heat,)

    ## Boundary condition
    Pin_test = 5.5      # inlet pressure (bar)
    yin_test = [1, 0]   # inlet composition (mol/mol)
    Tin_test = 300      # inlet temperature (K)
    Q_in_test = 0.02*0.031416*0.3  # volumetric flowrate (m^3/sec)
    Cvin_test = 5E-2    # inlet valve constant (m/sec/bar)
 
    Pout_test = 1       # outlet pressure (bar)
    Cvout_test = 2E-2   # outlet valve constant (m/sec/bar)
    c1.boundaryC_info(Pout_test,Pin_test,Tin_test,yin_test,
    Cvin_test,Cvout_test,Q_in_test,assigned_v_option = True)
 
    #c1.boundaryC_info(Pout_test, Pin_test,Tin_test,yin_test,Cvin_test)
 
    ## Initial condition
    P_init = 1*np.ones(N)                   # initial pressure (bar)
    y_init = [0*np.ones(N), 1*np.ones(N)]   # initial composition (mol/mol)
    T_init = 300*np.ones(N)                 # initial temperature (K)
    q_init = iso_fn_test([P_init*y_init[0],
    P_init*y_init[1]],T_init)               # initial uptake
    c1.initialC_info(P_init, T_init, T_init,y_init,q_init)
    
    ## print here
    print(c1)
    c1.run_mamo(10,CPUtime_print = True)
    c1.Graph(
        1,0,loc = [0.8,0.85], 
        yaxis_label=r'C$_{1}$ (mol/m$^{3}$)',
        figsize = [8.5,5])
    c1.Graph_P(
        1, loc = [1.1,0.85],
        figsize = [8.5,5],)
    plt.show( )
    