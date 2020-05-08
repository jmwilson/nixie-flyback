import math
import e_series

# Converter parameters
L = 10e-6 # H
V = 170. # V
V_g = 12. # V
n = 10
I_out = 30e-3 # A
f_s = 120e3 # Hz
I_pk = 3 # A, DA-2032 saturation current

R = V/I_out
D = (V/V_g)*math.sqrt(2*L*f_s/R)
D_CCM = V/(V + n*V_g)

# LM5155 controller parameters
V_CLTH = 100e-3 # V
V_SL = 40e-3 # V
I_SL = 30e-6 # A

# Derating factors
alpha = 1.2
beta = .6

# Compute R_s and R_SL
R_smax = 2*V_SL*f_s/(alpha * V/(n*L))
R_s0 = V_CLTH/I_pk

print(f"Without slope compensation: R_s = {R_s0:.3g}")
if R_s0 > R_smax:
    print("Slope compensation required in CCM")
    R_s = n*L*f_s*(D*V_SL + V_CLTH)/(D*V*beta + I_pk*n*L*f_s)
    R_SL = (beta*V*V_CLTH - I_pk*n*L*f_s*V_SL)/(I_SL*(beta*D*V + I_pk*n*L*f_s))
    print(f"Best choice: R_s = {R_s:.3g}, R_SL = {e_series.e24(R_SL)}")
    if R_SL > 2e3:
        print("R_SL exceeds maximum permitted by LM5155")
        R_SL = float(input("Choose R_SL between 0 and 2000: "))
        R_s = (V_CLTH - I_SL*R_SL*D)/I_pk
        print(f"Choice 1: meet current limit, R_s = {R_s:.3g}")
        print(f"Stability factor = {(V_SL + I_SL*R_SL)*f_s/(R_s*V/(n*L)):.3g}")

        R_s = (V_SL + I_SL*R_SL)*f_s/(beta*V/(n*L))
        print(f"Choice 2: meet stability requirement, R_s = {R_s:.3g}")
        print(f"DCM Current limit = {(V_CLTH - I_SL*R_SL*D)/R_s} A")
