import math
import numpy
import scipy
import e_series

# Converter parameters
L = 10e-6 # H
L_lk = 150e-9 # H
V = 170. # V
V_g = 5. # V
n = 10
I_out = 30e-3 # A
f_s = 350e3 # Hz

# Snubber parameters
V_sn = 32. # V
V_sn_ripple = 1. # V

# Peak current calculation
D = 1/(1 + n*V_g/V)
Dp = 1 - D
I_L = n/Dp*I_out
I_L_ripple = V_g/(2*L)*D/f_s
I_L_pk = I_L + I_L_ripple

# Snubber calculation
P = L_lk*I_L_pk**2/2*V_sn/(V_sn - V/n)*f_s
R_sn = V_sn**2/P
C_sn = V_sn/(V_sn_ripple*R_sn*f_s)

print("AN4147 RCD snubber")
print(f"Peak inductor current: {I_L_pk:.3g} A")
print(f"Power dissipated: {P:.3g} W")
print(f"R_sn = {e_series.e24(R_sn):.2g} Ω, C_sn = {e_series.e12(C_sn):.2g} F")
