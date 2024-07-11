import numpy as np
import math
from scipy.interpolate import splrep, splev, BSpline

class multiple_flux_surfaces():
    def __init__(self, curve_fitting_class_object):
        '''
        initites the needed variables 
        '''
        self.curve_fitting_class_object = curve_fitting_class_object
        
        # self.R_axis = R_axis
        # self.Z_axis = Z_axis
        # if (psi_fs_last is None):
        #     self.psi_fs_last = self.curve_fitting_class_object.flux_function.eval_psi(r_fs, z_fs)
        # else:
        #     self.psi_fs_last = self.curve_fitting_class_object.flux_function.psi_fs_last
        
        
        self.psi_curve_coefficients = None 
        self.R_ij = None
        self.Z_ij = None
        self.R_Bspline_triplets = None
        self.Z_Bspline_triplets = None
        self.dR_dTheta=None 
        self.dZ_dTheta = None
        self.dR_ds=None 
        self.dZ_ds = None
        self.rotational_transforms_psi = None
        self.Fourier_spline_triplets = None
        self.iota_spline_coeff = None
        
        self.R_spline_Fourier=None
        self.Z_spline_Fourier=None
    
    
    
    # def s_to_psi(self,s):
    #     '''
    #     converts from normalised radial coordinate s to psi value
    #     '''
    #     return (s)*(self.psi_fs_last-self.psi_axis)+self.psi_axis
    
    def get_multiple_flux_surfaces(self,N_s=20, N_theta=300, m_max_end=6):
        '''
        calculates the parametrisation corresponding to the flux surfaces at the normalised radial variable s with N_s equidistant points from 0 to 1
        '''
        S = np.linspace(0,1,N_s+1)
        R_ij = np.zeros((N_s+1,N_theta))
        Z_ij = np.zeros((N_s+1,N_theta))
        
        
        R_ij[0] = np.ones(N_theta)*self.curve_fitting_class_object.flux_function.R_axis
        Z_ij[0] = np.ones(N_theta)*self.curve_fitting_class_object.flux_function.Z_axis
        psi_curve_coefficients = [np.zeros(4*(m_max_end+1))]
        
        for i in range(1,N_s+1):
            #finding R value for each psi value
            # sol = opt.minimize(fun=self.psi_finder_dummy_function, args=psi_surfaces[i], x0=self.R_axis, tol=1.0e-15, # method='Nelder-Mead')
            
            #finding initial parametric curve
            # s_fs, x_out, xmap = self.curve_fitting_class_object.find_flux_surface(s_fs=S[i], #r_fs=sol.x[0], z_fs=self.Z_axis,   
            #                                                                         m_max_end=math.floor(10), m_max_start=2)
            
            #root finding using points from initial parametric curve
            # rr,zz= self.curve_fitting_class_object.get_curve(s_fs=s_fs,x_final_in=x_out,xmap_in=xmap,Np_in=2*(2*m_max_end+1)+1,rzcorr_in=1,thet_offset=0)
            
            # getting psi error for line search
            # del_s_points = np.average(np.abs(self.curve_fitting_class_object.flux_function.eval_s(rr,zz) - s_fs*np.ones(len(rr))))
            
            #inverse fourier transofr to get coefficients for final parametric curve
            # x_out,xmap = self.curve_fitting_class_object.inverse_fourier_tranform(m_max_end=m_max_end,rr=rr,zz=zz)
            
            x_out, xmap = self.curve_fitting_class_object.get_full_contour(s_fs=S[i], m_max_end=m_max_end)
            
            
            rr_curve,zz_curve=self.curve_fitting_class_object.eval_curve(xmap=xmap, x=x_out, N=1003)
            
            #psi error for final curve
            del_s_curve = np.average(np.abs(self.curve_fitting_class_object.flux_function.eval_s(rr_curve,zz_curve) - S[i]))
            
            
   
            
            
            #adding points and curve coefficients to respective lists/arrays
            psi_curve_coefficients.append((x_out,xmap))
            R_ij[i], Z_ij[i] = self.curve_fitting_class_object.eval_curve(xmap=xmap, x=x_out, N=N_theta)
            print(f"Final Curve (S={S[i]}) error: {del_s_curve} ")
        
        
        
        self.psi_curve_coefficients, self.R_ij, self.Z_ij = psi_curve_coefficients, R_ij, Z_ij
       
             
    
    
    def change_psi_fs_last(self, r_fs=None,z_fs=None,psi_fs_last=None, Recalculate=False):
        '''
        changes the last psi flux surface. If recalculate is true, it will recalculate all class atributes with the new last psi flux surface
        '''
        if (r_fs is not None) and (z_fs is not None):
            self.curve_fitting_class_object.flux_function.psi_fs_last = self.curve_fitting_class_object.flux_function.eval_psi(r_fs, z_fs)
        elif psi_fs_last is not None:
            self.curve_fitting_class_object.flux_function.psi_fs_last = psi_fs_last
            
        if Recalculate:
            self.get_multiple_flux_surfaces()
            self.get_BSplines()
            self.angular_derivative()
        
        
    def get_BSplines(self):
        '''
        gets coefficients for the radial Bsplies interpolation of the flux surfaces
        '''
    
        self.R_Bspline_triplets = []
        self.Z_Bspline_triplets = []
        
        for i in range(len(self.R_ij[0])):
            s= np.linspace(0,1,len(self.R_ij))
            self.R_Bspline_triplets.append(splrep(s,self.R_ij[:,i], k=3))
            self.Z_Bspline_triplets.append(splrep(s,self.Z_ij[:,i], k=3))
        
        self.R_Bspline_triplets = np.array(self.R_Bspline_triplets)
        self.Z_Bspline_triplets = np.array(self.Z_Bspline_triplets)
            
    
    
    
    
    
    
    def angular_derivative(self):
        '''
        calculates dr/dtheta, dz/dtheta of the flux surfaces using the fourier series parametrisation
        '''
        self.dR_dTheta = np.zeros(self.R_ij.shape)
        self.dZ_dTheta = np.zeros(self.Z_ij.shape)
        
        for i in range(1,self.R_ij.shape[0]):
            x0, xmap = self.psi_curve_coefficients[i]
            self.dR_dTheta[i], self.dZ_dTheta[i] = self.curve_fitting_class_object.eval_curve_dt(x=x0, xmap=xmap, N=self.R_ij.shape[1])
    
    
    def radial_derivative(self):
        '''
        calculates dr/dtheta, dz/dtheta of the flux surfaces using the BSplie interpolation
        '''
        self.dR_ds = np.zeros(self.R_ij.shape)
        self.dZ_ds = np.zeros(self.Z_ij.shape)
        s = np.linspace(0,1,self.R_ij.shape[0])
        for i in range(self.R_ij.shape[1]):
            
            self.dR_ds[:,i] = splev(x=s, tck=self.R_Bspline_triplets[i], der=1)
            self.dZ_ds[:,i] = splev(x=s, tck=self.Z_Bspline_triplets[i], der=1)
            
            
    
   
        
    
    
    
    
    def calculate_fourier_coefficients_splines_interpolate(self):
        all_coefficients = np.zeros((len(self.psi_curve_coefficients), len(self.psi_curve_coefficients[0])))
        xmap = 0
        
        
        for i in range(1,len(self.psi_curve_coefficients)):
            coeff, xmap = self.psi_curve_coefficients[i]
            all_coefficients[i] = coeff

        all_coefficients[0][xmap["str_r_c"]] = self.curve_fitting_class_object.flux_function.R_axis
        if(self.curve_fitting_class_object.unsymm):
            all_coefficients[0][xmap["str_z_c"]] = self.curve_fitting_class_object.flux_function.Z_axis  # be careful with unsymm=False...
        
        self.Fourier_spline_triplets = []
        
        s = np.linspace(0,1,len(self.psi_curve_coefficients))
        
        for i in range( len(self.psi_curve_coefficients[0])):
            self.Fourier_spline_triplets.append(splrep(s,all_coefficients[:,i], k=3))
    
    
    
    
    
    def sTheta_to_RZ(self, s, theta):
        
        xmap, coeffs = self.fourier_coefficients_spline_coefficients_unwrapper(s)
        
        return self.curve_fitting_class_object.eval_curve(x=coeffs, xmap=xmap, thet=theta)
    
    
    def angular_deriavtive_at_inputs(self,s,theta):
        xmap, coeffs = self.fourier_coefficients_spline_coefficients_unwrapper(s)
        
        return self.curve_fitting_class_object.eval_curve_dt(x=coeffs, xmap=xmap, thet=theta)
        
    
        
    
    def fourier_coefficients_spline_coefficients_unwrapper(self, s):
        coeffs = []
        
        for i in range(len(self.Fourier_spline_triplets)):
            t,c,k = self.Fourier_spline_triplets[i]
            coeffBSpline = BSpline(t,c,k, extrapolate=False)
            coeffs.append(coeffBSpline(s))
        
        
        coeffs = np.array(coeffs)
        xmap = self.psi_curve_coefficients[len(self.psi_curve_coefficients)-1][1]
        
        return xmap, coeffs
    
    
    def rotational_transform_using_mapping(self, s, N=200):
        #s=max(s,1.0e-4)
        
        F_s = self.curve_fitting_class_object.flux_function.eval_F(s)
        theta = np.linspace(0,1,N,endpoint=False)
        
        
        psi_axis = self.curve_fitting_class_object.flux_function.psi_axis
        psi_last = self.curve_fitting_class_object.flux_function.psi_fs_last
        
        R,Z = self.sTheta_to_RZ(s=s, theta=theta)
        dr_dtheta, dz_dtheta = self.angular_deriavtive_at_inputs(s=s, theta=theta)
        
        ds_dr, ds_dz = self.curve_fitting_class_object.flux_function.eval_del_s(R, Z)
        

        detf=(ds_dr*dz_dtheta-ds_dz*dr_dtheta)
        #assert all(detf>0.) ,  'detf must be >0 on the flux surface, curve parameter must be counter-clockwise!'+(("...min=%e,max=%e")%(np.amin(detf),np.amax(detf)))

        integrand = detf/(R * ((ds_dr)**2 + (ds_dz)**2))
        
        intermediate = (F_s/(2*np.pi*(psi_last - psi_axis))) * np.average(integrand)
        
        return 1./intermediate # <== this is iota
    
    
    def rotational_transform_with_dedicated_contour(self, s, N=400, **kwargs):
        # s=max(s,1.0e-4)
        if s <1.0e-4:
            return 1/self.safety_factor_taylor_approximation(s=s)
        
        F_s = self.curve_fitting_class_object.flux_function.eval_F(s)
        theta = np.linspace(0,1,N,endpoint=False)
        
        
        psi_axis = self.curve_fitting_class_object.flux_function.psi_axis
        psi_last = self.curve_fitting_class_object.flux_function.psi_fs_last
        
        curve_coeffs,curve_xmap= self.curve_fitting_class_object.get_full_contour(s_fs=s,**kwargs)
        
        R,Z = self.curve_fitting_class_object.eval_curve(x=curve_coeffs, xmap=curve_xmap,thet=theta)
        dr_dtheta, dz_dtheta = self.curve_fitting_class_object.eval_curve_dt(x=curve_coeffs, xmap=curve_xmap,thet=theta)
        
        ds_dr, ds_dz = self.curve_fitting_class_object.flux_function.eval_del_s(R, Z)
        

        detf=(ds_dr*dz_dtheta-ds_dz*dr_dtheta)
        
        assert all(detf>=0.) ,  'detf must be >0 on the flux surface, curve parameter must be counter-clockwise!'+(("...min=%e,max=%e")%(np.amin(detf),np.amax(detf)))

        integrand = detf/(R * ((ds_dr)**2 + (ds_dz)**2))
        
        intermediate = (F_s/(2*np.pi*(psi_last - psi_axis))) * np.average(integrand)
        
        return 1./intermediate  # <== this is iota
        
       
    
    def safety_factor_with_dedicated_contour(self, s, N=200, m_max_end = 12, **kwargs):
        return 1./(self.rotational_transform_with_dedicated_contour(s, N=N, m_max_end = m_max_end,**kwargs))
    
    
    def safety_factor_with_using_mapping(self, s, N=200):
        return 1./(self.rotational_transform_using_mapping(s, N=N))
    
    def safety_factor_taylor_approximation(self, s):
        Raxis, Zaxis = self.curve_fitting_class_object.flux_function.R_axis, self.curve_fitting_class_object.flux_function.Z_axis
        s_rr, s_zz = self.curve_fitting_class_object.flux_function.eval_del_squared_s(Raxis, Zaxis)
        psi_last = self.curve_fitting_class_object.flux_function.psi_fs_last
        psi_axis = self.curve_fitting_class_object.flux_function.psi_axis
        F_s = self.curve_fitting_class_object.flux_function.eval_F(s)
        
        # q_s = (F_s/(Raxis*(psi_last - psi_axis))) * (1/np.sqrt(s_rr*s_zz)) * (1 + s/(s_rr * Raxis**2) + (1.5)*(s/(s_rr * Raxis**2))**2)
        
        q_s = (F_s/(Raxis*(psi_last - psi_axis))) * (1/np.sqrt(s_rr*s_zz)) * (1/np.sqrt(1-2*s/(s_rr * (Raxis**2))))
        
        return q_s
    
    
    
    
    def rotational_transform_interpolate(self, N):
        
        s = np.linspace(0,1,N)

        iota = []

        for i in s:
            iota.append(self.rotational_transform_with_dedicated_contour(s=i))
            
        self.iota_spline_coeff = splrep(s,iota, k=3)
        
        