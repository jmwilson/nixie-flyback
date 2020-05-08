import math
import e_series

# Converter parameters
L = 10e-6 # H
V = 170. # V
V_g = 5. # V
n = 10
I_out = 30e-3 # A
f_s = 350e3 # Hz
I_pk = 3 # A, DA-2032 saturation current

D = V/(V + n*V_g)

# LM5155 controller parameters
V_CLTH = 100e-3 # V
V_SL = 40e-3 # V
I_SL = 30e-6 # A

# Derating factors
alpha = 1.2
beta = .82

# Compute R_s and R_SL
R_smax = 2*V_SL*f_s/(alpha * V/(n*L))
R_s0 = V_CLTH/I_pk

print(f"Without slope compensation: R_s = {R_s0:.3g}, maximum permissible R_s = {R_smax:.3g}")
if R_s0 > R_smax:
	R_s = n*L*f_s*(D*V_SL + V_CLTH)/(D*V*beta + I_pk*n*L*f_s)
	R_SL = (beta*V*V_CLTH - I_pk*n*L*f_s*V_SL)/(I_SL*(beta*D*V + I_pk*n*L*f_s))
	print("Slope compensation required")
	print(f"R_s = {R_s:.3g}, R_SL = {e_series.e24(R_SL)}")
