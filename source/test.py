from scipy.interpolate import splev
import numpy as np

x = splev(x=np.linspace(0,1,20), tck=(3,,3))