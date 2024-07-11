from curve_fitting.curve_fitting_class import curve_fitter
from matplotlib import pyplot as plt
from grad_shafranov_class import GS_Solution
import numpy as np

fit1 = curve_fitter(formfactor="Dshape", solverClass=GS_Solution)

# plotting norm vector(gradient of psi)
"""Rgrid, Zgrid = fit1.Profile_solver.get_RZ_grid(20)
B_field_R, B_field_Z = fit1.Profile_solver.eval_del_psi(Rgrid, Zgrid)
plt.quiver(Rgrid, Zgrid, B_field_R, B_field_Z)
plt.show()"""

# plotting psi contourmap
r_fs = 3
psi_fs = fit1.Profile_solver.eval_psi(r_fs, 0)
Rgrid, Zgrid = fit1.Profile_solver.get_RZ_grid(150)
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(Rgrid, Zgrid, fit1.Profile_solver.eval_psi(Rgrid, Zgrid), levels=40, cmap='jet', alpha=0.9)

r0 = fit1.Profile_solver.mid_poit_finder()
ax.plot(r0[0], r0[1], 'ro')

cp2 = ax.contour(Rgrid, Zgrid, fit1.Profile_solver.eval_psi(Rgrid, Zgrid), [psi_fs], colors='k', linewidths=2.,
                 linestyles='solid')


fig.colorbar(cp)  # Add a colorbar to a plot
ax.set_title('Psi functions evaluated(R,Z)')
ax.set_xlabel('R')
ax.set_ylabel('Z')
plt.show()
# root searching and the plotting results
psi_fs, results = fit1.find_flux_surface(r_fs=r_fs, m_max_start=2, m_max_end=8)
res=results["m_max=6"]
x_final,xmap=res["x_final"],res["xmap"]
rr,zz= fit1.get_curve(psi_fs=psi_fs,x_final_in=x_final,xmap_in=xmap,Np_in=1024,rzcorr_in=1,thet_offset=0)
plt.plot(rr, zz, 'bo')
print("psi fs: ", psi_fs)
del_psi = np.average(np.abs(fit1.Profile_solver.eval_psi(rr,zz) - psi_fs*np.ones(len(rr))))
print("del_psi: " , del_psi)
plt.xlabel('R')
ax.axis("equal")
plt.show()


