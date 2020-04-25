import math
import e_series

# Converter parameters
L = 10e-6 # H
V = 170. # V
V_g = 5. # V
n = 10
I_out = 30e-3 # A
f_s = 350e3 # Hz

D = V/(V + n*V_g)
Dp = 1 - D
I_L = n/Dp * I_out
I_Lripple = V_g/(2*L)*D/f_s
I_pk = I_L + I_Lripple
print(f"I_L = {I_L:.3g}, I_pk = {I_pk:.3g}")

# LM5155 controller parameters
V_CLTH = 100e-3 # V
V_SL = 40e-3 # V
I_SL = 30e-6 # A

# Derating factors
alpha = .2
beta = .82

# Compute R_s and R_SL
R_smax = (1 - alpha)*2*V_SL*f_s*n*L/V
R_s0 = (V_CLTH - V_SL*D)/I_pk

print(f"Without slope compensation: R_s = {R_s0:.3g}, maximum permissible R_s = {R_smax:.3g}")
if R_s0 >= R_smax:
	R_s = n*L*V_CLTH*f_s/(D*V*beta + (1 + alpha)*I_pk*n*L*f_s)
	R_SL = (beta*V*(V_CLTH - D*V_SL) - (1 + alpha)*n*L*V_SL*I_pk*f_s)/(I_SL*(beta*D*V + (1 + alpha)*I_pk*n*L*f_s))
	print("Slope compensation required")
	print(f"R_s = {R_s:.3g}, R_SL = {e_series.e24(R_SL)}")
