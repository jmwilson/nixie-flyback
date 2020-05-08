import math
import numpy
import scipy
import e_series
from ltiarithmetic import Z_R, Z_C, Z_L, TransferFunction
import matplotlib.pyplot as plt

def transfer_plot(system):
    freq = numpy.logspace(1,6,30 * 6)
    plt.style.use('dark_background')
    fig, ax = plt.subplots(constrained_layout=True)
    fig.set_size_inches(10,6.18)
    fig.set_dpi(150)
    tau = 2 * math.pi
    w, mag, phase = system.bode(tau * freq)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("gain (dB)")
    h1, = ax.semilogx(w/tau, mag, label="gain")
    ax.grid(True, which="both", axis="x", color="xkcd:dark grey")
    phase_ax = ax.twinx()
    phase_ax.set_ylabel("phase (deg)")
    h2, = phase_ax.semilogx(w/tau, phase, linestyle="dashed", label="phase")
    plt.legend(handles=[h1, h2], loc="upper right")
    return fig

# LM5155 error amplifier characteristics
g_m = 2e-3  # S
G_comp = 0.142

# LM5155 current mode characteristics
V_sl = 40e-3 # V
I_sl = 30e-6 # A
V_fb = 1. # V

# Design parameters
V = 170. # V
V_g = 9 # V
f_s = 120e3 # Hz
L = 10e-6 # H
n = 10.
R_f = 33e-3 # Ohm
C = 2*220e-9 # F
R_sl = 0. # Ohm
I_out = 30e-3 # A

# Desired crossover frequency and phase margin
f_c = 12e3 # Hz
p_m = 60. # degrees

R = V/I_out # Ohm

# Current mode slopes (A/s)
M_a = (V_sl + I_sl*R_sl)*f_s/R_f
M_1 = V_g/L

# Find best resistive divider
max_error = math.inf
for r1 in e_series.e96_values:
    r1 *= 10**6
    r2 = e_series.e96(r1 / (V - V_fb))
    error = abs(r2/(r1 + r2) - V_fb/V)
    if error < max_error:
        max_error = error
        R_1, R_2 = r1, r2

# G_vc: control to output gain
Gc0 = math.sqrt(2*L*f_s/R)*1/(1 + M_a/M_1)*R/2
wp = 2/(R*C)
G_vc = Gc0 * TransferFunction([1], [1/wp, 1])

# H: resistive divider sensor gain
H = R_2 / (R_1 + R_2)

# T_u: Uncompensated loop gain
T_u = H * G_comp * 1/R_f * G_vc

# Find uncompensated gain at desired crossover and get starting point for R_c
_, [T_u_f_c] = T_u.freqresp([2 * math.pi * f_c])
gain_f_c_u = numpy.abs(T_u_f_c)
R_c_0 = 1/(g_m * gain_f_c_u)

# Get starting point for C_c1
phase_f_c_u = numpy.angle(T_u_f_c, deg=True)
required_phase_adj = (p_m - 180.) - phase_f_c_u
f_z = f_c * math.exp(math.pi/2 * required_phase_adj/45)
C_c1_0 = 1/(2 * math.pi * f_z * R_c_0)

def par(*args):
    return 1/sum(map(lambda x: 1/x, args))

# Optimize for gain=1 and desired phase margin at crossover frequency f_c
def loop_system(x, f, margin):
    r_c, c_c1 = x
    G_c = g_m * (Z_R(r_c) + Z_C(c_c1))
    system = H * G_c * G_comp * 1/R_f * G_vc
    _, [resp] = system.freqresp([2 * math.pi * f])
    gain = numpy.abs(resp)
    phase_margin = 180 + numpy.angle(resp, deg=True)
    return (gain - 1, phase_margin - margin)

[R_c, C_c1] = scipy.optimize.fsolve(
    loop_system,
    [R_c_0, C_c1_0],
    (f_c, p_m)
)

print("DCM flyback type II OTA compensator")

print(f"Initial estimate: R_c_0 = {R_c_0:.2g} Ω, C_c1_0 = {C_c1_0:.2g} F")
print(
    f"Compensation network for crossover at {f_c:g} Hz with {p_m:g}°"
    f" phase margin: R_c = {R_c:g} Ω, C_c1 = {C_c1:g} F"
)

R_c = e_series.e24(R_c)
C_c1 = e_series.e6(C_c1)

# G_c: compensated error amplifier gain
G_c = g_m * (Z_R(R_c) + Z_C(C_c1))

# T: compensated loop gain
T = H * G_c * G_comp * 1/R_f * G_vc

def system_crossover(f, system):
    _, [resp] = system.freqresp([2 * math.pi * f])
    return numpy.abs(resp) - 1

[crossover] = scipy.optimize.fsolve(system_crossover, f_c, T)
_, [T_crossover] = T.freqresp([2 * math.pi * crossover])
phase_margin = numpy.angle(T_crossover, deg=True) + 180.
print(f"Component values: R_c = {R_c:.2g} Ω, C_c1 = {C_c1:.2g} F")
print(f"Crossover at {crossover:.2g} Hz, phase margin: {phase_margin:.2g}°")

transfer_plot(G_c).suptitle("Compensator frequency response")
transfer_plot(T).suptitle("Loop frequency response")
plt.show()
