# convergence checker
# for rk3 scheme w/ 32 32 10 (nx ny nz) we got the order of 2.9999524917331124
import numpy as np
qx_coarse = float(8.5266663293100242e-1)
qx_mid = float(8.5235636071483245e-1)
qx_fine =  float(8.5231757541062358e-1)


p = np.log((qx_mid-qx_coarse)/(qx_fine-qx_mid))/np.log(2.0)
print(p)

qx_coarse = float(8.6014596418435896e-1)
qx_mid = float(8.4627259677896149e-1)
qx_fine =  float(8.2979463661100616e-1)

p = np.log((qx_mid-qx_coarse)/(qx_fine-qx_mid))/np.log(2.0)
print(p)
# 2nd order flat wall  w/ rk3  
qx_coarse = float( 4.8068146233492687E-01)
qx_mid = float(4.6808135220443342E-01)
qx_fine =  float(4.6891703492944486E-01)

p_2nd = np.log((qx_mid-qx_coarse)/(qx_fine-qx_mid))/np.log(2.0)
print(p_2nd)