<<<<<<< HEAD
from grad_shafranov_class import GS_Solution
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from scipy import integrate


class curve_fitter:
    def __init__(self, R0=0, p_eps=0, p_kappa=0, p_delta=0, p_A=0, qaxis=0, p0=0, formfactor=None,
                 solverClass=GS_Solution):
        if formfactor is None:
            self.Profile_solver = solverClass(R0=R0, p_eps=p_eps, p_kappa=p_kappa, p_delta=p_delta, p_A=p_A,
                                              qaxis=qaxis,
                                              p0=p0)
        else:
            self.Profile_solver = solverClass(formfactor=formfactor)

        self.parametric_curve = 0
=======
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


class curve_fitter:
    def __init__(self, solverClassObject,unsymm=True):
        self.flux_function = solverClassObject
        self.unsymm=unsymm

        
>>>>>>> iota_without_s_derivative

    def eval_1d_fourier(self, t, coef_c=[], coef_s=[]):
        '''Evaluates cosine/sine series, first coef_c is mode m=0, of coef_s is mode m=1'''
        x = 0.0 * t
        for m in range(0, len(coef_c)):
            x += coef_c[m] * np.cos(m * t * 2 * np.pi)
        for m in range(0, len(coef_s)):
<<<<<<< HEAD
            x += coef_s[m] * np.sin((m + 1) * t * 2 * np.pi)
        return x
=======
            x += coef_s[m] * np.sin((m) * t * 2 * np.pi)
        return x
    
    def eval_1d_fourier_coeff_gradient(self, t, m_c=[], m_s=[], R = True):
        
        
        T_c,M_c = np.meshgrid(t,m_c)
        T_s,M_s = np.meshgrid(t,m_s)
        
        grad = np.concatenate((np.cos(M_c * T_c * 2 * np.pi),np.sin(M_s * T_s *2 * np.pi)))
        add_on = np.zeros(grad.shape)
        if R:
            return np.concatenate((grad,add_on))
        else:
            return np.concatenate((add_on,grad))
        
        
        
>>>>>>> iota_without_s_derivative

    def eval_1d_fourier_dt(self, t, coef_c=[], coef_s=[]):
        '''Evaluates first derivative of cosine/sine series, first coef_c is mode m=0, of coef_s is mode m=1'''
        dxdt = 0.0 * t
        for m in range(0, len(coef_c)):
            dxdt += -m * 2 * np.pi * coef_c[m] * np.sin(m * t * 2 * np.pi)
        for m in range(0, len(coef_s)):
<<<<<<< HEAD
            dxdt += (m + 1) * 2 * np.pi * coef_s[m] * np.cos((m + 1) * t * 2 * np.pi)

        return dxdt

    def eval_curve(self, x, xmap, N, thet_offset=0.):
=======
            dxdt += (m) * 2 * np.pi * coef_s[m] * np.cos((m) * t * 2 * np.pi)

        return dxdt
    
    
    def eval_1d_fourier_dt_coeff_gradient(self, t, m_c=[], m_s=[], R=True):
        T_c,M_c = np.meshgrid(t,m_c)
        T_s,M_s = np.meshgrid(t,m_s)
        
        grad = np.concatenate((-M_c * 2 * np.pi*np.sin(M_c * T_c * 2 * np.pi),M_s * 2 * np.pi * np.cos(M_s * T_s *2 * np.pi)))
        add_on = np.zeros(grad.shape)
        if R:
            return np.concatenate((grad,add_on))
        else:
            return np.concatenate((add_on,grad))
    
    

    def eval_curve(self, x, xmap, N=100, thet_offset=0., thet=None):
>>>>>>> iota_without_s_derivative
        '''
        Evaluates r,z of the curve defined by fourier coefficients
        rcoef_c/s and zcoef_c/s , at N points. coefficients are mapped from 1d array x, using xmap dict.
        '''
<<<<<<< HEAD
        thet = np.linspace(0. + thet_offset, 1. + thet_offset, N, endpoint=False)
=======
        if thet is None:
            thet = np.linspace(0. + thet_offset, 1. + thet_offset, N, endpoint=False)
        
>>>>>>> iota_without_s_derivative
        # thet=thet+0.5*(thet[1]-thet[0]) #shift by 1/2 dx
        rr = self.eval_1d_fourier(thet,
                                  coef_c=x[xmap["str_r_c"]:xmap["str_r_c"] + xmap["nr_c"]],
                                  coef_s=x[xmap["str_r_s"]:xmap["str_r_s"] + xmap["nr_s"]])

        zz = self.eval_1d_fourier(thet,
                                  coef_c=x[xmap["str_z_c"]:xmap["str_z_c"] + xmap["nz_c"]],
                                  coef_s=x[xmap["str_z_s"]:xmap["str_z_s"] + xmap["nz_s"]])

        return rr, zz
<<<<<<< HEAD

    def eval_curve_dt(self, x, xmap, N, thet_offset=0.):
=======
    
    
    def eval_curve_coefficient_gradient(self, x, xmap, N=100, thet_offset=0., thet=None):
        
        x = xmap["m"]
        if thet is None:
            thet = np.linspace(0. + thet_offset, 1. + thet_offset, N, endpoint=False)
        
        grad_rr = self.eval_1d_fourier_coeff_gradient(thet,
                                  m_c=x[xmap["str_r_c"]:xmap["str_r_c"] + xmap["nr_c"]],
                                  m_s=x[xmap["str_r_s"]:xmap["str_r_s"] + xmap["nr_s"]], R=True)
        
        grad_zz = self.eval_1d_fourier_coeff_gradient(thet,
                                  m_c=x[xmap["str_z_c"]:xmap["str_z_c"] + xmap["nz_c"]],
                                  m_s=x[xmap["str_z_s"]:xmap["str_z_s"] + xmap["nz_s"]], R=False)
        
        return grad_rr, grad_zz
        

    def eval_curve_dt(self, x, xmap, N=100, thet_offset=0., thet=None):
>>>>>>> iota_without_s_derivative
        '''
        Evaluates r,z of the curve defined by fourier coefficients
        rcoef_c/s and zcoef_c/s , at N points. coefficients are mapped from 1d array x, using xmap dict.
        '''
<<<<<<< HEAD
        thet = np.linspace(0. + thet_offset, 1. + thet_offset, N, endpoint=False)
=======
        
        if thet is None:
            thet = np.linspace(0. + thet_offset, 1. + thet_offset, N, endpoint=False)
        
>>>>>>> iota_without_s_derivative
        # thet=thet+0.5*(thet[1]-thet[0]) #shift by 1/2 dx
        dr_dt = self.eval_1d_fourier_dt(thet,
                                        coef_c=x[xmap["str_r_c"]:xmap["str_r_c"] + xmap["nr_c"]],
                                        coef_s=x[xmap["str_r_s"]:xmap["str_r_s"] + xmap["nr_s"]])

        dz_dt = self.eval_1d_fourier_dt(thet,
                                        coef_c=x[xmap["str_z_c"]:xmap["str_z_c"] + xmap["nz_c"]],
                                        coef_s=x[xmap["str_z_s"]:xmap["str_z_s"] + xmap["nz_s"]])
        return dr_dt, dz_dt
<<<<<<< HEAD

    def eval_curve_normal(self, x, xmap, N, **kwargs):
        dr_dt, dz_dt = self.eval_curve_dt(x, xmap, N, **kwargs)
        nn = np.sqrt(dr_dt ** 2 + dz_dt ** 2)
        return dz_dt / nn, -dr_dt / nn  # n_r,n_z

    def psidiff(self, x_in, xmap, psi_goal, N):
        '''
        function for the minimizer: evaluate psi on the curve and
        return the squared normalized difference to psi_goal
        '''
        rr, zz = self.eval_curve(x_in, xmap, N)
        psi_curve = self.Profile_solver.eval_psi(rr, zz)
        return (psi_curve - psi_goal) ** 2 / (psi_goal ** 2)

    def BdotN(self, x_in, xmap, N):
        # weird format ask pls
        rr, zz = self.eval_curve(x_in, xmap, N)
        Bpol_r, Bpol_z = self.Profile_solver.eval_Bpol(rr, zz)
=======
    
    
    def eval_curve_dt_coefficient_gradient(self, x, xmap, N=100, thet_offset=0., thet=None):
        
        x = xmap["m"]
        
        if thet is None:
            thet = np.linspace(0. + thet_offset, 1. + thet_offset, N, endpoint=False)
        
        grad_dr_dt = self.eval_1d_fourier_dt_coeff_gradient(thet,
                                  m_c=x[xmap["str_r_c"]:xmap["str_r_c"] + xmap["nr_c"]],
                                  m_s=x[xmap["str_r_s"]:xmap["str_r_s"] + xmap["nr_s"]], R=True)
        
        grad_dz_dt = self.eval_1d_fourier_dt_coeff_gradient(thet,
                                  m_c=x[xmap["str_z_c"]:xmap["str_z_c"] + xmap["nz_c"]],
                                  m_s=x[xmap["str_z_s"]:xmap["str_z_s"] + xmap["nz_s"]], R=False)
        
        return grad_dr_dt, grad_dz_dt
    

    def eval_curve_normal(self, x, xmap, N,**kwargs):
        '''
        returns normal of the curve
        '''
        dr_dt, dz_dt = self.eval_curve_dt(x, xmap, N, **kwargs)
        
        
        
        nn = np.sqrt(dr_dt ** 2 + dz_dt ** 2)
        return dz_dt / nn, -dr_dt / nn  # n_r,n_z

    def S_diff(self, x_in, xmap, s_goal, N):
        '''
        function for the minimizer: evaluate s on the curve and
        return the squared normalized difference to s_goal
        '''
        rr, zz = self.eval_curve(x_in, xmap, N)
        s_curve = self.flux_function.eval_s(rr, zz)
        return (s_curve - s_goal) ** 2 
    
    def S_diff_coeff_gradient(self, x_in, xmap, s_goal, N):
        rr, zz = self.eval_curve(x_in, xmap, N)
        s_curve = self.flux_function.eval_s(rr, zz)
        
        grad_rr, grad_zz = self.eval_curve_coefficient_gradient(x_in, xmap, N)
         
        ds_dr, ds_dz = self.flux_function.eval_del_s(rr,zz)
        
        summation_array = 2 * (s_curve - s_goal) * (ds_dr*grad_rr + ds_dz * grad_zz)
        
        return summation_array
    

    def BdotN(self, x_in, xmap, N):
        '''
        returns dot product of the B field vector and the vector normal to the parametrised curve
        '''
        rr, zz = self.eval_curve(x_in, xmap, N)
        Bpol_r, Bpol_z = self.flux_function.eval_Bpol(rr, zz)
>>>>>>> iota_without_s_derivative

        dr_dt, dz_dt = self.eval_curve_dt(x_in, xmap, N)

        n_r = -dz_dt
        n_z = dr_dt
        return (n_r * Bpol_r + n_z * Bpol_z) ** 2 / ((Bpol_r ** 2 + Bpol_z ** 2) * (n_r ** 2 + n_z ** 2))
<<<<<<< HEAD

    def eval_curve_normal(self, x, xmap, N, **kwargs):
        dr_dt, dz_dt = self.eval_curve_dt(x, xmap, N, **kwargs)
        nn = np.sqrt(dr_dt ** 2 + dz_dt ** 2)
        return dz_dt / nn, -dr_dt / nn  # n_r,n_z
=======
    
    
    def BdotN_coeff_gradient(self, x_in, xmap, N):
        rr, zz = self.eval_curve(x_in, xmap, N)
        Bpol_r, Bpol_z = self.flux_function.eval_Bpol(rr, zz)

        dr_dt, dz_dt = self.eval_curve_dt(x_in, xmap, N)
        R0 = self.flux_function.R0
        
        grad_dr_dt, grad_dz_dt = self.eval_curve_dt_coefficient_gradient(x_in, xmap, N)
        
        grad_rr ,grad_zz = self.eval_curve_coefficient_gradient(x_in, xmap, N)
        
        d2psi_dR2, d2psi_dZ2, d2psi_dRdZ = self.flux_function.eval_del_squared_psi(R=rr,Z=zz,Mixed=True)
        
        v = ((Bpol_r ** 2 + Bpol_z ** 2) * (dz_dt ** 2 + dr_dt ** 2))
        u = (-dz_dt * Bpol_r + dr_dt * Bpol_z) ** 2
        
        grad_Bpol_r = (-1/R0)*(d2psi_dRdZ*grad_rr + d2psi_dZ2*grad_zz)
        grad_Bpol_z = (1/R0)*(d2psi_dR2*grad_rr + d2psi_dRdZ*grad_zz)
        
        du = 2*(-dz_dt * Bpol_r + dr_dt * Bpol_z)*(-(dz_dt*grad_Bpol_r + grad_dz_dt*Bpol_r) + (dr_dt*grad_Bpol_z + grad_dr_dt*Bpol_z))
        dv = 2*((Bpol_r ** 2 + Bpol_z ** 2)*(dz_dt*grad_dz_dt + dr_dt*grad_dz_dt) + (dz_dt ** 2 + dr_dt ** 2)*(Bpol_r*grad_Bpol_r + Bpol_z*grad_Bpol_z))
        
        return (v*du - u*dv)/(v**2)
    


    def grads_dot_t(self, x_in, xmap, N):
        '''
        dot product of the tangential of the curve with the normal of the flux surface, normalized to the cosine^2.
        returns , (ds/dR*dR/dt + ds/dZ*dZ/dt)^2/((ds/dR^2+ds/dZ^2)*(dR/dt^2+dZ/dt^2)
        '''
        rr, zz = self.eval_curve(x_in, xmap, N)

        ds_dR,ds_dZ = self.flux_function.eval_del_s(rr, zz)

        dR_dt, dZ_dt = self.eval_curve_dt(x_in, xmap, N)

        return (ds_dR * dR_dt + ds_dZ * dZ_dt) ** 2 / (( ds_dR** 2 + ds_dZ ** 2) * (dR_dt ** 2 + dZ_dt ** 2))

>>>>>>> iota_without_s_derivative

    def Mscale(self, x_in, xmap, p=1, q=1):

        return np.sum(xmap["m"] ** (p + q) * x_in ** 2) / np.sum(xmap["m"] ** (p) * x_in ** 2)
<<<<<<< HEAD

    def minf(self, x_in, xmap, psi_goal, N):
        '''
        function for the minimizer:
        '''
        # return 0.5*np.sum(psidiff(x_in,xmap,psi_goal,N))/N
        return 0.5 * np.sum(self.psidiff(x_in, xmap, psi_goal, N)) / N + 0.5 * np.sum(
            self.BdotN(x_in, xmap, N)) / N  # <==== add normal condition
=======
    
    def Mscale_grad(self, x_in, xmap, p=1, q=1):
        u =np.sum(xmap["m"] ** (p + q) * x_in ** 2)
        v= np.sum(xmap["m"] ** (p) * x_in ** 2)
        
        return 2*x_in*(v*(xmap["m"]**(p+q)) - u*(xmap["m"]**(p)))/(v**2)

    def minf(self, x_in, xmap, s_goal, N, Mscale_factor=1.0e-11, p=2):
        '''
        function for the minimizer:
        ''' 
        
        # return 0.5*np.sum(s_diff(x_in,xmap,s_goal,N))/N
        return (1-Mscale_factor)*(np.sum(self.S_diff(x_in, xmap, s_goal, N)) / N + np.sum(
            self.BdotN(x_in, xmap, N)) / N  )+ Mscale_factor*self.Mscale(x_in, xmap,p=p) # <==== add normal condition + Mscale
        
    
    def minf_gradient(self, x_in, xmap, s_goal, N, Mscale_factor=1.0e-11, p=2):
        '''
        function for gradient with respect to curve coefficients(x_in) of minimizer function:]
        ''' 
        
        
        # print(f"S grad:{(np.sum(self.S_diff_coeff_gradient(x_in, xmap, s_goal, N),1))/N}\n, Mscale grad:{Mscale_factor*self.Mscale_grad(x_in, xmap,p=p)}\n, Bdot Grad:{ np.sum(self.BdotN_coeff_gradient(x_in, xmap, N),1)/N}")
        
        
        return (1-Mscale_factor)*(np.sum(self.S_diff_coeff_gradient(x_in, xmap, s_goal, N),1))/N + Mscale_factor*self.Mscale_grad(x_in, xmap,p=p) + np.sum(self.BdotN_coeff_gradient(x_in, xmap, N),1)/N
    
    
>>>>>>> iota_without_s_derivative

    # def minf_LS(x_in,xmap,psi_goal,N):
    #    return psidiff(x_in,xmap,psi_goal,N) + BdotN(x_in,xmap,N)

<<<<<<< HEAD
    def minf_point(self, x, psi_goal):
        return 0.5 * (self.Profile_solver.eval_psi(x[0], x[1]) - psi_goal) ** 2 / (psi_goal ** 2)
=======
    def minf_point(self, x, s_goal):
        return 0.5 * (self.flux_function.eval_s(x[0], x[1]) - s_goal) ** 2 
>>>>>>> iota_without_s_derivative

    def eval_line(self, a, xp, xs):
        return [xp[0] + a * xs[0], xp[1] + a * xs[1]]

    def minf_line(self, a, xp, xs, psi_goal):
        return self.minf_point(self.eval_line(a, xp, xs), psi_goal)

<<<<<<< HEAD
    def psi_line(self, a, xp, xs, psi_goal):
        x = self.eval_line(a, xp, xs)
        return self.Profile_solver.eval_psi(x[0], x[1]) - psi_goal
=======
    def S_line(self, a, xp, xs, s_goal):
        x = self.eval_line(a, xp, xs)
        return self.flux_function.eval_s(x[0], x[1]) - s_goal
>>>>>>> iota_without_s_derivative

    def point_distance(self, x_in, xmap, N, rr_ref, zz_ref):
        rr, zz = self.eval_curve(x_in, xmap, N)
        return (rr - rr_ref) ** 2 + (zz - zz_ref) ** 2

    def minf_dist(self, x_in, xmap, N, rr_ref, zz_ref):
        return 0.5 * np.sum(self.point_distance(x_in, xmap, N, rr_ref, zz_ref)) / N
<<<<<<< HEAD

    def find_flux_surface(self, r_fs=3.15, z_fs=0., m_max_start=4, m_max_end=6, **kwargs):
=======
    
    
    def find_bounding_box(self, s_fs):
        '''
        Finds box that bounds the particular psi contour
        '''
        R, Z = self.flux_function.get_RZ_grid(gridPoints=150)
        splot = self.flux_function.eval_s(R, Z)
        
        s_bb =max(s_fs,0.01)   # safe side bounding box close to magnetic axis
 
        fig, ax = plt.subplots(figsize=(10, 10))

        im = ax.contourf(R, Z, splot, 40, cmap='jet', alpha=0.4)
        ax.contour(R, Z, splot, 40, cmap='jet')
        
        ax.contour(R, Z, splot, [s_bb], colors='k', linewidths=2., linestyles='solid')

        # get the contour object
        fs = ax.contour(R, Z, splot, [s_bb])

        fs_list = fs.collections[0].get_paths()
        # Sort the paths by its length. Assume main one is the longest(?)
        fs_list.sort(key=len, reverse=True)

        fs_coord = fs_list[0].vertices
        fs_r = fs_coord[:, 0]
        fs_z = fs_coord[:, 1]

        idx_left = np.min(fs_r[:])
        idx_right = np.max(fs_r[:])
        idx_lower = np.min(fs_z[:])
        idx_upper = np.max(fs_z[:])

        # print ("bound. box, left (%f,%f),right (%f,%f),lower (%f,%f),upper (%f,%f)" % (fs_r[idx_left],fs_z[idx_left],fs_r[idx_right],fs_z[idx_right],fs_r[idx_lower],fs_z[idx_lower],fs_r[idx_upper],fs_z[idx_upper]))
        # ax.plot(fs_r[0:-1:10],fs_z[0:-1:10],'bx')
        
        plt.close()
        
        return idx_left, idx_right, idx_lower, idx_upper
        

    def find_flux_surface(self, r_fs=None, z_fs=None, s_fs=None, m_max_start=4, m_max_end=6, **kwargs):
>>>>>>> iota_without_s_derivative
        '''
        Main function to find a closed flux surface (with parametrization of the angle!) of the coil+plasma field.
        Visualization helps if default parameters would change.
        Defines a contour(=flux surface) of psi at a given point (r_fs,z_fs), uses the contour object and extract the closed contour.
        The contour object gives a bounding box which is used to initialize the fourier series of the curve.
        Then the coefficients are optimized to match the contour level,
        by minimizing the difference of the flux evaluated at the curve to the flux of the contour.
        '''
        from scipy.optimize import minimize
        from scipy.optimize import root
<<<<<<< HEAD

        R, Z = self.Profile_solver.get_RZ_grid(gridPoints=150)
        psi = self.Profile_solver.eval_psi(R, Z)

        # select flux surface going through the point:

        psi_fs = self.Profile_solver.eval_psi(r_fs, z_fs)

        fig, ax = plt.subplots(figsize=(10, 10))

        im = ax.contourf(R, Z, psi, 40, cmap='jet', alpha=0.4)
        ax.contour(R, Z, psi, 40, cmap='jet')

        ax.contour(R, Z, psi, [psi_fs], colors='k', linewidths=2., linestyles='solid')

        # get the contour object
        fs = ax.contour(R, Z, psi, [psi_fs])

        fs_list = fs.collections[0].get_paths()
        # Sort the paths by its length. Assume main one is the longest(?)
        fs_list.sort(key=len, reverse=True)

        fs_coord = fs_list[0].vertices
        fs_r = fs_coord[:, 0]
        fs_z = fs_coord[:, 1]

        idx_left = np.argmin(fs_r[:])
        idx_right = np.argmax(fs_r[:])
        idx_lower = np.argmin(fs_z[:])
        idx_upper = np.argmax(fs_z[:])

        # print ("bound. box, left (%f,%f),right (%f,%f),lower (%f,%f),upper (%f,%f)" % (fs_r[idx_left],fs_z[idx_left],fs_r[idx_right],fs_z[idx_right],fs_r[idx_lower],fs_z[idx_lower],fs_r[idx_upper],fs_z[idx_upper]))
        # ax.plot(fs_r[0:-1:10],fs_z[0:-1:10],'bx')

        psi_back = self.Profile_solver.eval_psi(fs_r, fs_z)

        print("psi_fs %e, max abs (psi_contour-psi_fs)= %e" % (psi_fs, np.amax(np.abs(psi_back - psi_fs))))

        # testing minimizer with one point:
        x0_point = [3.5, 1.2]
        diff_psi = self.minf_point(x0_point, psi_fs)
        print("psi_fs %e, initial residual,sqrt||(psi_point-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))

        res_point = minimize(self.minf_point, x0_point, args=(psi_fs), tol=1.0e-16)

        print("message minimizer_point: " + res_point.message)
        diff_psi = self.minf_point(res_point.x, psi_fs)
        print("psi_fs %e, final residual,sqrt||(psi_point-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))
        print("final diff_psi |psi_point-psi_fs|= %e " % (np.sqrt(2 * psi_fs ** 2 * diff_psi)))
        # plt.plot([x0_point[0],res_point.x[0]],[x0_point[1],res_point.x[1]],'ro')

        # testing minimizer with point along line:
        x0_a = 0.1;
        xp = [3.5, 1.2];
        xs = [1, -0.5]
        diff_psi = self.minf_line(x0_a, xp, xs, psi_fs)
        print("psi_fs %e, initial residual,sqrt||(psi_line-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))
        res_line = minimize(self.minf_line, x0_a, args=(xp, xs, psi_fs), tol=1.0e-16)
        print("message minimizer_line: " + res_line.message)
        diff_psi = self.minf_line(res_line.x, xp, xs, psi_fs)
        print("psi_fs %e, final residual,sqrt||(psi_line-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))
=======
        from scipy.optimize import LinearConstraint

        # select flux surface going through the point:
        if(s_fs is None):
          assert((not(r_fs is None)) and (not (z_fs is None)))
          s_fs = self.flux_function.eval_s(r_fs, z_fs)
             

        idx_left,idx_right,idx_lower,idx_upper = self.find_bounding_box(s_fs=s_fs)

        #print("psi_fs %e, max abs (psi_contour-psi_fs)= %e" % (psi_fs, np.amax(np.abs(psi_back - psi_fs))))

        # testing minimizer with one point:
        # x0_point = [3.5, 1.2]
        # diff_psi = self.minf_point(x0_point, psi_fs)
        # print("psi_fs %e, initial residual,sqrt||(psi_point-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))

        # res_point = minimize(self.minf_point, x0_point, args=(psi_fs), tol=1.0e-16)

        # print("message minimizer_point: " + res_point.message)
        # diff_psi = self.minf_point(res_point.x, psi_fs)
        # print("psi_fs %e, final residual,sqrt||(psi_point-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))
        # print("final diff_psi |psi_point-psi_fs|= %e " % (np.sqrt(2 * psi_fs ** 2 * diff_psi)))
        # # plt.plot([x0_point[0],res_point.x[0]],[x0_point[1],res_point.x[1]],'ro')

        # # testing minimizer with point along line:
        # x0_a = 0.1
        # xp = [3.5, 1.2]
        # xs = [1, -0.5]
        # diff_psi = self.minf_line(x0_a, xp, xs, psi_fs)
        # print("psi_fs %e, initial residual,sqrt||(psi_line-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))
        # res_line = minimize(self.minf_line, x0_a, args=(xp, xs, psi_fs), tol=1.0e-16)
        # print("message minimizer_line: " + res_line.message)
        # diff_psi = self.minf_line(res_line.x, xp, xs, psi_fs)
        # print("psi_fs %e, final residual,sqrt||(psi_line-psi_fs)^2/psi_fs^2||= %e " % (psi_fs, np.sqrt(diff_psi)))
>>>>>>> iota_without_s_derivative

        #################
        # find the contour parametrization
        #################
<<<<<<< HEAD
        all_results = {}
        unsymm = True
        for m_max in range(m_max_start, m_max_end + 1):
            print("=====m_max= %d ====== " % (m_max))
            rcoef_c0 = np.zeros(m_max + 1)
            if (unsymm):
                rcoef_s0 = np.zeros(m_max)  # unsymmetric
=======
        #all_results = {}
        xmap=0
        x_out=0
        
        for m_max in range(m_max_start, m_max_end + 1):
            #print("=====m_max= %d ====== " % (m_max))
            rcoef_c0 = np.zeros(m_max + 1)
            if (self.unsymm):
                rcoef_s0 = np.zeros(m_max+1)  # unsymmetric
>>>>>>> iota_without_s_derivative
                zcoef_c0 = np.zeros(m_max + 1)  # unsymmetric
            else:
                rcoef_s0 = []  # up-down symmetric
                zcoef_c0 = []  # up-down symmetric
<<<<<<< HEAD
            zcoef_s0 = np.zeros(m_max)
=======
            zcoef_s0 = np.zeros(m_max+1)
>>>>>>> iota_without_s_derivative

            # first initialization:
            # initialize  up-down symmetric curve inside bounding box of contour
            if (m_max == m_max_start):
                # rcoef_c0[2]=0.1 #m=2
<<<<<<< HEAD
                rcoef_c0[1] = 0.4 * (fs_r[idx_right] - fs_r[idx_left])  # m=1
                rcoef_c0[0] = 0.5 * (fs_r[idx_right] + fs_r[idx_left])  # - rcoef_c0[2] # m=0
                zcoef_s0[0] = 0.4 * (fs_z[idx_upper] - fs_z[idx_lower])  # m=1
                zcoef_c0[0] = 0.5 * (fs_z[idx_upper] + fs_z[idx_lower])
                # zcoef_s0[1]=-0.1 #m=2
            # use previous solution
            else:
                rcoef_c0[0:m_max] = res.x[xmap["str_r_c"]:xmap["str_r_c"] + xmap["nr_c"]]
                if (unsymm):
                    rcoef_s0[0:m_max - 1] = res.x[xmap["str_r_s"]:xmap["str_r_s"] + xmap["nr_s"]]
                    zcoef_c0[0:m_max] = res.x[xmap["str_z_c"]:xmap["str_z_c"] + xmap["nz_c"]]
                zcoef_s0[0:m_max - 1] = res.x[xmap["str_z_s"]:xmap["str_z_s"] + xmap["nz_s"]]

            # BUILD 1D solution vector x0
            x0 = np.concatenate((rcoef_c0, rcoef_s0, zcoef_c0, zcoef_s0))
            xmap = {}
            xmap["nr_c"] = len(rcoef_c0)
            xmap["nr_s"] = len(rcoef_s0)
            xmap["nz_c"] = len(zcoef_c0)
            xmap["nz_s"] = len(zcoef_s0)
            xmap["str_r_c"] = 0
            xmap["str_r_s"] = 0 + xmap["nr_c"]
            xmap["str_z_c"] = xmap["str_r_s"] + xmap["nr_s"]
            xmap["str_z_s"] = xmap["str_z_c"] + xmap["nz_c"]
            # mode number
            xmap["m"] = np.concatenate((np.arange(0, xmap["nr_c"]),
                                        np.arange(0, xmap["nr_s"]) + 1,
                                        np.arange(0, xmap["nz_c"]),
                                        np.arange(0, xmap["nz_s"]) + 1))
=======
                rcoef_c0[1] = 0.4 * (idx_right - idx_left)  # m=1
                rcoef_c0[0] = 0.5 * (idx_right + idx_left)  # - rcoef_c0[2] # m=0
                zcoef_s0[1] = 0.4 * (idx_upper - idx_lower)  # m=1
                if(self.unsymm):
                  zcoef_c0[0] = 0.5 * (idx_upper + idx_lower)
                # zcoef_s0[1]=-0.1 #m=2
            # use previous solution
            else:
                
                if (self.unsymm):
                    rcoef_s0[0:m_max] = res.x[xmap["str_r_s"]:xmap["str_r_s"] + xmap["nr_s"]]
                    zcoef_c0[0:m_max] = res.x[xmap["str_z_c"]:xmap["str_z_c"] + xmap["nz_c"]]
                zcoef_s0[0:m_max] = res.x[xmap["str_z_s"]:xmap["str_z_s"] + xmap["nz_s"]]
                rcoef_c0[0:m_max] = res.x[xmap["str_r_c"]:xmap["str_r_c"] + xmap["nr_c"]]

            # BUILD 1D solution vector x0
            x0 = np.concatenate((rcoef_c0, rcoef_s0, zcoef_c0, zcoef_s0))
            xmap = self.build_xmap(len(rcoef_c0), len(rcoef_s0), len(zcoef_c0), len(zcoef_s0))
>>>>>>> iota_without_s_derivative

            # number of points for evaluating the "distance" in psi on the curve
            N = 1 + 4 * m_max  # add another +1 if theta starts at 1/2*dx
            N_post = 251

            rr, zz = self.eval_curve(x0, xmap, N)

            # plt.plot(rr,zz,'r.')

<<<<<<< HEAD
            diff_psi = self.minf(x0, xmap, psi_fs, N_post)
            print("psi_fs %e, initial residual= %e " % (psi_fs, diff_psi))

            res = minimize(self.minf, x0, args=(xmap, psi_fs, N), tol=1.0e-10)
            print("message minimizer: " + res.message)
=======
            diff_s = self.minf(x0, xmap, s_fs, N_post)
            #print("psi_fs %e, initial residual= %e " % (psi_fs, diff_psi))
            A=np.zeros((2,len(x0)))
            A[0,xmap["str_r_s"]]=1.
            A[1,xmap["str_z_s"]]=1.
            constraint = LinearConstraint(A,lb=0.,ub=0.)
            "constraints=constraint,"
            res = minimize(self.minf, x0, args=(xmap, s_fs, N), tol=1.0e-14, constraints=constraint, jac=self.minf_gradient)
            #print("message minimizer: " + res.message)
>>>>>>> iota_without_s_derivative
            x_out = res.x

            # res_LS=least_squares(minf_LS, x0, args=(xmap,psi_fs,N))
            # print("message minimizer LS: " + res_LS.message)
            # x_out=res_LS.x

            # correct points onto contour and do another least squares fit:
            # rr_ref,zz_ref=eval_curve(x_out,xmap,N)
            # nr_ref,nz_ref=eval_curve_normal(x_out,xmap,N)
            # for i in range(0,N):
            #    x0_a=0.
            #    xp=[rr_ref[i],zz_ref[i]]
            #    xs=[nr_ref[i],nz_ref[i]]
            #    res_line=root(psi_line, x0_a, args=(xp,xs,psi_fs),tol=1.0e-14)
            #    rr_ref[i],zz_ref[i]=eval_line(res_line.x,xp,xs)
            # psi_back=eval_full_psi(rr_ref,zz_ref)
            # print ("fit points: min/max (psi_fit-psi_fs)= %e %e" % (np.amin(psi_back-psi_fs),np.amax(psi_back-psi_fs)))

            # res_dist=minimize(minf_dist, x_out, args=(xmap,N,rr_ref,zz_ref),tol=1.0e-14)
            # print("message minimizer distance: " + res_dist.message)
            # x_out=res_dist.x

            # post-processing

<<<<<<< HEAD
            diff_psi = self.minf(x_out, xmap, psi_fs, N_post)
            print("psi_fs %e, final residual = %e " % (psi_fs, diff_psi))

            rr, zz = self.eval_curve(x_out, xmap, N_post)
            psi_back = self.Profile_solver.eval_psi(rr, zz)
            print("psi_fs %e, min/max (psi_fit-psi_fs)= %e %e" % (
                psi_fs, np.amin(psi_back - psi_fs), np.amax(psi_back - psi_fs)))

            print("maximum error in |B.n|/|B|= %e" % (np.amax(np.abs(np.sqrt(self.BdotN(x_out, xmap, N_post))))))
            print("Mscale(p=1/2/4)= %f %f %f" % (
                self.Mscale(x_out, xmap, p=1), self.Mscale(x_out, xmap, p=2), self.Mscale(x_out, xmap, p=4)))

            # add result
            which_m = ("m_max=%d" % (xmap["nz_s"]))
            all_results[which_m] = {"x_final": x_out, "xmap": xmap}
=======
            diff_s = self.minf(x_out, xmap, s_fs, N_post)
            #print("s_fs %e, final residual = %e " % (s_fs, diff_s))

            rr, zz = self.eval_curve(x_out, xmap, N_post)
            s_back = self.flux_function.eval_s(rr, zz)
            # print("s_fs %e, min/max (s_fit-s_fs)= %e %e" % (
            #     s_fs, np.amin(s_back - s_fs), np.amax(s_back - s_fs)))

            # print("maximum error in |B.n|/|B|= %e" % (np.amax(np.abs(np.sqrt(self.BdotN(x_out, xmap, N_post))))))
            # print("Mscale(p=1/2/4)= %f %f %f" % (
            #     self.Mscale(x_out, xmap, p=1), self.Mscale(x_out, xmap, p=2), self.Mscale(x_out, xmap, p=4)))

            # add result
            #which_m = ("m_max=%d" % (xmap["nz_s"]))
            #all_results[which_m] = {"x_final": x_out, "xmap": xmap}
>>>>>>> iota_without_s_derivative
        # visualize result
        # # rr, zz = self.eval_curve(x_out, xmap, 50)
        # plt.plot(rr, zz, 'bo')

        # # visualize normalized Bpol
        # # Bpol_r,Bpol_z = eval_full_Bpol(rr,zz)
        # # absBpol=np.sqrt(Bpol_r**2+Bpol_z**2)
        # # ax.quiver(rr,zz,Bpol_r,Bpol_z)

        # plt.xlabel('R')
        # #plt.ylabel('Z')
        # # plt.title('Plot of the poloidal field with "plasma" ')
        # plt.title((
        #         'Plot of the reference poloidal field (coils+ "plasma") \n with fitted closed flux surface through point (%4.2f,%4.2f)' % (
        #     r_fs, z_fs)))
        # fig.colorbar(im)
        # ax.axis("equal")
        # plt.show()
<<<<<<< HEAD

        return psi_fs, all_results

    def get_curve(self, psi_fs, x_final_in, xmap_in, Np_in, rzcorr_in, thet_offset=0):
=======
        # print("psi_fs %e, min/max (psi_fit-psi_fs)= %e %e" % (
        #       psi_fs, np.amin(psi_back - psi_fs), np.amax(psi_back - psi_fs)))

        # print("maximum error in |B.n|/|B|= %e" % (np.amax(np.abs(np.sqrt(self.BdotN(x_out, xmap, N_post))))))
        # print("Mscale(p=1/2/4)= %f %f %f" % (
        # self.Mscale(x_out, xmap, p=1), self.Mscale(x_out, xmap, p=2), self.Mscale(x_out, xmap, p=4)))
        
        return s_fs, x_out, xmap

    def get_curve(self, s_fs, x_final_in, xmap_in, Np_in, rzcorr_in, thet_offset=0):
        '''
        takes points from the curve, and perfofmrs a root search along the line normal to the curve.
        This is to incrase the accuary of the parametrisation
        '''
>>>>>>> iota_without_s_derivative
        from scipy.optimize import root
        rr, zz = self.eval_curve(x_final_in, xmap_in, Np_in, thet_offset=thet_offset)
        if rzcorr_in == 1:
            # set points exactly onto the psi_fs contour, via linesearch from curve point into normal direction
            rr_ref, zz_ref = self.eval_curve(x_final_in, xmap_in, Np_in, thet_offset=thet_offset)
            nr_ref, nz_ref = self.eval_curve_normal(x_final_in, xmap_in, Np_in, thet_offset=thet_offset)
            for i in range(0, Np_in):
                x0_a = 0.
                xp = [rr_ref[i], zz_ref[i]]
                xs = [nr_ref[i], nz_ref[i]]
<<<<<<< HEAD
                res_line = root(self.psi_line, x0_a, args=(xp, xs, psi_fs), tol=1.0e-19)
=======
                res_line = root(self.S_line, x0_a, args=(xp, xs, s_fs), tol=1.0e-19)
>>>>>>> iota_without_s_derivative
                rr[i], zz[i] = self.eval_line(res_line.x, xp, xs)
        return rr, zz
    
    
<<<<<<< HEAD
    def final_fourier_coefficients(self, N_points, m_max_end, rr, zz):
        r_c,r_s,z_c,z_s = [],[],[],[]
        r0,z0 = self.fourier_integrate(rr,zz,0)
        r_c.append(r0)
        z_c.append(z0)
        # for i in range():
        
         
    
    
    def fourier_integrate(self,rr,zz,m):
        exp = lambda M, t: np.exp(complex(0,1)*t*M*2*np.pi)
        
        x = np.linspace(0,1,len(rr), endpoint=False)
        RR = exp(m,x)*rr
        ZZ = exp(m,x) * zz
        
        Integral_rr = integrate.cumtrapz(x=x, y=RR, initial=0)
        Integral_zz = integrate.cumtrapz(x=x, y=ZZ, initial=0)
        
        
        if m == 0:
            cr_m = Integral_rr
            cz_m = Integral_zz
            return np.real(cr_m), np.real(cz_m)
        else:
            cr_m = 2*Integral_rr
            cz_m = 2*Integral_zz
            return np.real(cr_m), np.imag(cr_m), np.real(cz_m), np.imag(cz_m)
=======
    def inverse_fourier_tranform(self,m_max_end, rr, zz):
        '''
        performes inverse fourier transform of a given array of (rr,zz)points 
        '''
        rcoef_c0 , rcoef_s0 = self.inverse_fourier_tranform_1d(m_max_end=m_max_end, input=rr)
        zcoef_c0 , zcoef_s0 = self.inverse_fourier_tranform_1d(m_max_end=m_max_end, input=zz) 
            
        x0 = np.concatenate((rcoef_c0, rcoef_s0, zcoef_c0, zcoef_s0))
        xmap = self.build_xmap(len(rcoef_c0), len(rcoef_s0), len(zcoef_c0), len(zcoef_s0))
        
        return x0, xmap
    
    
    def get_full_contour(self, r_fs=None, z_fs=None, s_fs=None, m_max_start_first=1, m_max_end_first=6,rzcorr_in=1,thet_offset=0, m_max_end=12,root_search =True, **kwargs):
        assert m_max_end > 0
        psi_fs,lcfs_x,xmap=self.find_flux_surface(s_fs=s_fs, m_max_start=m_max_start_first, m_max_end=m_max_end_first)
        if not root_search:
            return lcfs_x,xmap
        
        rr,zz= self.get_curve(s_fs=s_fs,x_final_in=lcfs_x,xmap_in=xmap,Np_in=2*(2*m_max_end+1)+1,rzcorr_in=rzcorr_in,thet_offset=thet_offset)
        # print(f"N_in{2*(2*m_max_end+1)+1}")
        
        x_corr,xmap_corr = self.inverse_fourier_tranform(m_max_end=m_max_end,rr=rr,zz=zz)
        
        return x_corr,xmap_corr
      
        
    def inverse_fourier_tranform_1d(self,m_max_end, input):
        
        coef_s0 = np.zeros(m_max_end+1)
        coef_c0 = np.zeros(m_max_end+1)
        for i in range(m_max_end+1):
            coef_c0[i], coef_s0[i]   = self.fourier_integrate_1d(input=input,m=i)
        
        return coef_c0,  coef_s0
    
    
    def fourier_integrate_1d(self,input,m):
        '''
        integrates over an array to find the fourier coefficient corresponding to a specified harmonic m
        '''
        exp = lambda M, t: np.exp(complex(0,1)*t*M*2*np.pi)
        
        x = np.linspace(0,1,len(input), endpoint=False)
        RR = exp(m,x)*input
        
        #trapezoidal rule on periodic equidistant data on interval [0,1]
        Integral_rr = np.average(RR)
        
        
        if m == 0:
            return np.real(Integral_rr),0
        else:
            cr_m = 2*Integral_rr
            return np.real(cr_m), np.imag(cr_m)
        
        
    def build_xmap(self, lrc,lrs,lzc,lzs):
        '''
        builds xmap dictionary used for evaluating the curve from the fourier coefficients
        '''
        xmap = {}
        xmap["nr_c"] = lrc
        xmap["nr_s"] = lrs
        xmap["nz_c"] = lzc
        xmap["nz_s"] = lzs
        xmap["str_r_c"] = 0
        xmap["str_r_s"] = 0 + xmap["nr_c"]
        xmap["str_z_c"] = xmap["str_r_s"] + xmap["nr_s"]
        xmap["str_z_s"] = xmap["str_z_c"] + xmap["nz_c"]
        # mode number
        xmap["m"] = np.concatenate((np.arange(0, xmap["nr_c"]),
                                    np.arange(0, xmap["nr_s"]),
                                    np.arange(0, xmap["nz_c"]),
                                    np.arange(0, xmap["nz_s"])))
        
        return xmap
            
>>>>>>> iota_without_s_derivative
            
            
        
