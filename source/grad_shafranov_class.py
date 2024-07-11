from sympy import *
import numpy as np
import scipy.optimize as opt


class GS_Solution:
    def __init__(self, R0=0, p_eps=0, p_kappa=0, p_delta=0, p_A=0, qaxis=0, p0=0, formfactor=None, B0=1., psi_fs_last =0.):
        if formfactor is None:
            self.R0, self.p_eps, self.p_kappa, self.p_delta, self.p_A, self.qaxis, self.p0 = [R0, p_eps, p_kappa,
                                                                                              p_delta, p_A, qaxis, p0]
        else:
            dictionary = {
                "circular shaped tokamak": [5.0, 0.2, 1.0, 0.0, 0.0, 1.6, 0.02],
                "ITER from Cerfon": [6.2, 0.32, 1.7, 0.33, -0.155, 1.6, 0.08],
                "spherical tokamak NSTX": [1, 0.78, 2, 0.35, 0.0, 1.1, 0.001],
                "Dshape": [5.0, 0.32, 1.7, 0.33, -0.1, 1.6, 0.08]
            }
            self.R0, self.p_eps, self.p_kappa, self.p_delta, self.p_A, self.qaxis, self.p0 = dictionary[formfactor]
        
        self.psi_fs_last = psi_fs_last
        self.psi_profile = 0
        self.psi_function = 0
        self.DelPsiExpression = 0
        self.DelPsiFunc = 0
        self.Del2PsiFunc = 0
        self.B0 = B0
        self.F0 =0.
        self.dF2ds =0.
        self.psi_0 =0.
        self.R_axis = 0.
        self.Z_axis = 0.
        self.psi_axis=0.
        self.dpds =0.
        self.flux_function_solver()
        
        
    def s_to_psi(self,s):
        '''
        converts from normalised radial coordinate s to psi value
        '''
        return (s)*(self.psi_fs_last-self.psi_axis) + self.psi_axis

    def flux_function_solver(self):
        '''
        solving for psi flux functon
        '''
        # defining sympy symbols
        x = symbols('x', real=True,positive=True)
        y = symbols('y', real=True)
        eps, kappa, delta = symbols('epsilon, kappa, delta', real=True, positive=True)
        A = symbols('A', real=True)
        c1, c2, c3, c4, c5, c6, c7 = symbols('c1, c2, c3, c4, c5, c6, c7', real=True)
        alpha = asin(delta)

        # Defining basis functions
        psi0 = x ** 4 / 8 + A * (1 / 2 * x ** 2 * ln(x) - x ** 4 / 8)
        psi1 = 1.
        psi2 = x ** 2
        psi3 = y ** 2 - x ** 2 * ln(x)
        psi4 = x ** 4 - 4 * x ** 2 * y ** 2
        psi5 = 2 * y ** 4 - 9 * y ** 2 * x ** 2 + 3 * x ** 4 * ln(x) - 12 * x ** 2 * y ** 2 * ln(x)
        psi6 = x ** 6 - 12 * x ** 4 * y ** 2 + 8 * x ** 2 * y ** 4
        psi7 = 8 * y ** 6 - 140 * y ** 4 * x ** 2 + 75 * y ** 2 * x ** 4 - 15 * x ** 6 * ln(x) \
               + 180 * x ** 4 * y ** 2 * ln(x) - 120 * x ** 2 * y ** 4 * ln(x)

        # horner(psi1.diff(y,2).subs(ln(x),ln(xx)))
        xx = symbols('xx', real=True)
        tmp0 = horner(psi0.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))
        tmp1 = 0.  # horner(psi1.diff(y,2).subs(ln(x),ln(xx)))
        tmp2 = horner(psi2.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))
        tmp3 = horner(psi3.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))
        tmp4 = horner(psi4.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))
        tmp5 = horner(psi5.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))
        tmp6 = horner(psi6.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))
        tmp7 = horner(psi7.diff(x, 1).diff(y, 1).subs(ln(x), ln(xx)))

        psix = psi0 + c1 * psi1 + c2 * psi2 + c3 * psi3 + c4 * psi4 + c5 * psi5 + c6 * psi6 + c7 * psi7
        # psi

        # N coefficients for Boundary constraints
        N1 = - (1 + alpha) ** 2 / (eps * kappa ** 2)
        N2 = + (1 - alpha) ** 2 / (eps * kappa ** 2)
        N3 = - kappa / (eps * cos(alpha) ** 2)

        # setting up boundary constraint equations
        eqs = [Eq(S(psix.subs(x, 1 + eps).subs(y, 0)), S(0)), Eq(S(psix.subs(x, 1 - eps).subs(y, 0)), S(0)),
               Eq(S(psix.subs(x, 1 - delta * eps).subs(y, kappa * eps)), S(0)),
               Eq(S(psix.diff(x).subs(x, 1 - delta * eps).subs(y, kappa * eps)), S(0)),
               Eq(S(psix.diff(y, 2).subs(x, 1 + eps).subs(y, 0) + N1 * psix.diff(x).subs(x, 1 + eps).subs(y, 0)), S(0)),
               Eq(S(psix.diff(y, 2).subs(x, 1 - eps).subs(y, 0) + N2 * psix.diff(x).subs(x, 1 - eps).subs(y, 0)), S(0)),
               Eq(S(psix.diff(x, 2).subs(x, 1 - delta * eps).subs(y, kappa * eps)) + N3 * psix.diff(y).subs(x,
                                                                                                            1 - delta * eps).subs(
                   y, kappa * eps), S(0))]

        eqsx = []

        for eq in eqs:
            eqsx.append(simplify(
                eq.subs(eps, self.p_eps).subs(kappa, self.p_kappa).subs(delta, self.p_delta).subs(A, self.p_A)))

        constants = solve(eqsx, [c1, c2, c3, c4, c5, c6, c7])

        psi = simplify(
            psix.subs(eps, self.p_eps).subs(kappa, self.p_kappa).subs(delta, self.p_delta).subs(A, self.p_A).subs(c1,
                                                                                                                  constants[
                                                                                                                      c1]).subs(
                c2, constants[c2]).subs(c3, constants[c3]).subs(c4, constants[c4]).subs(c5, constants[c5]).subs(c6,
                                                                                                                constants[
                                                                                                                    c6]).subs(
                c7, constants[c7]))

        self.psi_profile = psi
        self.psi_function = lambdify((x, y), self.psi_profile, "numpy")
        self.DelPsiExpression = np.array([self.psi_profile.diff(x), self.psi_profile.diff(y)])
        self.DelPsiFunc = [lambdify((x, y), self.psi_profile.diff(x), "numpy"),
                           lambdify((x, y), self.psi_profile.diff(y), "numpy")]
        
        self.Del2PsiFunc =[lambdify((x, y), self.psi_profile.diff(x).diff(x), "numpy"),
                           lambdify((x, y), self.psi_profile.diff(y).diff(y), "numpy"),
                           lambdify((x, y), self.psi_profile.diff(x).diff(y), "numpy")]
        
        self.defining_important_values()
       
    
    
    def defining_important_values(self):
         # q_axis=F_axis/(R_axis*(d_rr(psireal)*d_zz(psireal))^1/2), psireal=psi*psi0, evaluated at axis...
        #       =F_axis/psi0 * 1/(R_axis*(d_rr(psi)*d_zz(psi))^1/2), from book Freidberg, ideal MHD, p.134
        # R=x*R0,Z=y*R0
        # d/dR=(d/dx)*1/R0, d/dZ=(d/dy)*1/R0
        xaxis, yaxis = self.find_axis_xy(1,0)
        
        dpsi_dx, dpsi_dy = self.DelPsiFunc
        d2psi_dr2, d2psi_dz2, dpsi_dxdy = self.Del2PsiFunc
        
        print( f"dpsi_dx:{dpsi_dx(xaxis, yaxis)}, dpsi_dy:{dpsi_dy(xaxis, yaxis)}")
        assert (np.allclose([dpsi_dx(xaxis, yaxis), dpsi_dy(xaxis, yaxis)], [0,0])) , f'gradient at axis NOT ZERO!!!'
        
        print(f'dpsi_dxdy:{dpsi_dxdy(xaxis, yaxis)}')
        assert (np.allclose(dpsi_dxdy(xaxis, yaxis), 0)), f'mixed derivative at axis NOT ZERO!!!'
        
        # print(f"xaxis= {xaxis},  yaxis= {yaxis},  non-nomralised psi = {self.psi_function(xaxis, yaxis)}")
        
       
        
        self.R_axis,self.Z_axis = xaxis*self.R0, yaxis*self.R0
        
        # print(f"R_axis= {self.R_axis},  Z_axis= {self.Z_axis}")
        
        self.F0=self.B0*self.R0
       
        # print(f"F0= {self.F0}")
        
        dpsi_dxdx, dpsi_dydy, dpsi_dxdy = self.Del2PsiFunc
        

        
        
        psi_dxdx_axis=dpsi_dxdx(xaxis,yaxis)
        psi_dydy_axis=dpsi_dydy(xaxis,yaxis)
        
        # print(f"D^2psi/dx^2= {psi_dxdx_axis}  D^2psi_dy^2= {psi_dydy_axis}")
        
        qaxis_fact=self.R0/((xaxis)*(psi_dxdx_axis*psi_dydy_axis)**0.5) #  R0*xaxis*(diff^2/R0^2*diff^2/R0^2)^1/2, 
        
        # print(f"qaxis fact= {qaxis_fact}")
        
        
        
        self.psi_0=qaxis_fact*self.F0/self.qaxis
        
        # print(f"psi0= {self.psi_0}")
        
        
        self.psi_axis = self.eval_psi(self.R_axis,self.Z_axis)
       
       
        # print(f"psi_axis= {self.psi_axis}")
        

        #B0=F0/R0  => F0=B0*R0
        
        #from qaxis, F0~Faxis, psi0=F_axis/q_axis * (1/R_axis...) = R0/q_axis * (1/R_axis...) 
        # F^2=(F_0^2+dF^2/ds *s)  is linear in s~psi, from Cerfon&Freidberg paper 
        #  d/dpsi(F^2) = 2*F*d/dpsi(F) =-2*A*psi0^2/R0^2, 
        #  s=(psi-psi_axis)/(psi_edge-psi_axis), psi_edge=0 #normalized radial coordinate
        #  ds/dpsi = -1/psi_axis
        # -> d/ds(F^2)=-(2*A*psi0^2/R0^2) *(psi_last-psi_axis) ,
        self.dF2ds = -2*self.p_A*self.psi_0**2*(self.psi_fs_last-self.psi_axis)/(self.R0**2)
        
        #  print(f"d^2F/ds: {self.dF2ds}")
        
        #gradient of the pressure over normalized flux, p=p0+dpds*s 
        #mu0=4*pi*1.0e-07 
        mu0=1
        self.dpds=(1-self.p_A)*self.psi_0**2*self.psi_axis/(mu0*self.R0**4)
        
    
    def eval_F(self,s):
        '''
        Inputs: psi (scalar or numpy array)
        
        Output: Free Function F evaluated at the inputted psi value/ array of psi values
        '''
        return np.sqrt(self.F0**2+self.dF2ds*s)
    
    def eval_p(self,s):
        '''
        Inputs: 
        psi (scalar or numpy array)
        
        Output: 
        P function values evaluated at the inputted psi value/ array of psi values
        '''
        return (self.p0+self.dpds*s)

    def eval_s(self, R, Z):
        ''' s~psi, s=0 at axis and s=1 at the boundary '''
        return (self.eval_psi(R, Z)-self.psi_axis)/(self.psi_fs_last-self.psi_axis)

    def eval_psi(self, R, Z):
        '''
        Evaluared psi function
        
        Inputs:
        R: 2d numpy array
        Z: 2d numpy array
        
        Outputs
        psi function values evaluated at R,Z: 2d aray
        '''
        X, Y = [R / self.R0, Z / self.R0]
        return self.psi_0*self.psi_function(X, Y)

    def get_RZ_grid(self, gridPoints, Rlim=0, Zlim=0, dd_scaling=1.2):
        '''
        Returns Mesh grid for R and Z
        
        Inputs:
        gridPoints(int)
        Rlim(float)
        Zlim(float)
        
        Outputs:
        Rgrid(2d numpy array)
        Zgrid(2d numpy array)
        '''
        if Rlim == 0 and Zlim == 0:
            dd = 1.7 * self.p_eps * self.p_kappa * self.R0
            R = np.linspace(max(self.R0 - dd, 10 ** -2), self.R0 + dd, gridPoints)
            Z = np.linspace(-1.3 * dd, +1.3 * dd, gridPoints)
            Rgrid, Zgrid = np.meshgrid(R, Z)
            return Rgrid, Zgrid
        else:
            Nr = gridPoints
            Nz = int((Zlim[1] - Zlim[0]) / (Rlim[1]-Rlim[0]))*Nr
            r = np.linspace(Rlim[0], Rlim[1], Nr)
            z = np.linspace(Zlim[0], Zlim[1], Nz)
            R, Z = np.meshgrid(r, z)
            return R, Z

    def eval_del_psi(self, R, Z):
        '''
        evaluates gradient of psi function
        
        Inputs:
        R(2d numpy array)
        Z(2d numpy array)
        
        Outputs:
        dpsi_dR(2d numpy array)
        dpsi_dZ(2d numpy array)
        '''
        dpsi_dx, dpsi_dy = self.DelPsiFunc
        X, Y = [R / self.R0, Z / self.R0]
        dpsi_dR, dpsi_dZ = [dpsi_dx(X, Y) *self.psi_0/ self.R0, dpsi_dy(X, Y) *self.psi_0/ self.R0]
        return dpsi_dR, dpsi_dZ
    
    
    def eval_del_squared_psi(self,R,Z, Mixed=False):
        '''
        evaluates d^2/dR^2 and d^2/dZ^2 of psi function
        
        Inputs:
        R(2d numpy array)
        Z(2d numpy array)
        
        Outputs:
        d2psi_dR2(2d numpy array)
        d2psi_dZ2 2d numpy array)
        '''
        
        d2psi_dxx, d2psi_dyy, d2psi_dxy = self.Del2PsiFunc
        X, Y = [R / self.R0, Z / self.R0]
        d2psi_dR2, d2psi_dZ2, d2psi_dRdZ = [d2psi_dxx(X, Y) *self.psi_0/(self.R0**2), d2psi_dyy(X, Y) *self.psi_0/ (self.R0**2), d2psi_dxy(X, Y) *self.psi_0/ (self.R0**2)]
        if Mixed==True:
            return d2psi_dR2, d2psi_dZ2, d2psi_dRdZ
        
        return d2psi_dR2, d2psi_dZ2
    
    
    
    def eval_del_s(self, R, Z):
        '''
        evaluates gradient of psi function
        
        Inputs:
        R(2d numpy array)
        Z(2d numpy array)
        
        Outputs:
        ds_dR(numpy array)
        ds_dZ(numpy array)
        '''
        return (self.eval_del_psi(R=R, Z=Z) / (self.psi_fs_last-self.psi_axis))
    
    
    def eval_del_squared_s(self,R,Z):
        '''
        evaluates d^2/dR^2 and d^2/dZ^2 of s function
        
        Inputs:
        R(2d numpy array)
        Z(2d numpy array)
        
        Outputs:
        d2s_dR2(2d numpy array)
        d2s_dZ2 2d numpy array)
        '''
        
        return (self.eval_del_squared_psi(R=R, Z=Z) / ((self.psi_fs_last-self.psi_axis)))
        
        
    def find_axis_xy(self,xinit, yinit):
        sol = opt.minimize(fun=self.psi_xy_dummy_func, x0=[xinit,yinit], tol=1.0e-10, jac=self.psi_xy_dummy_gradient_vector)
        return sol.x

    def psi_xy_dummy_func(self, X):
        ''' Function to minimize in find_axis_xy'''
        return self.psi_function(X[0], X[1])
    
    def psi_xy_dummy_gradient_vector(self, X):
        ''' gradient in xy to minimize in find_axis_xy'''
        dpsi_dx, dpsi_dy = self.DelPsiFunc
        return [dpsi_dx(X[0],X[1]) , dpsi_dy(X[0],X[1])]


    def eval_Bpol(self, R, Z):
        '''
        Evaluates B_pol
        
        Inputs:
        R(2d numpy array)
        Z(2d numpy array)
        
        Outputs:
        B_pol(1d list)
        '''
        dpsi_dr, dpsi_dz = self.eval_del_psi(R, Z)
        return -dpsi_dz / self.R0, dpsi_dr / self.R0
    
    
    
    
    
