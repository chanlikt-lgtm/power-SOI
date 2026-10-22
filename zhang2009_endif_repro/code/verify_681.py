"""Verify the paper's claim  "E_I = 681 V/um for sigma_in = 4e12, t_S=1um, eps_I=3.9eps0".

Three independent checks + a hypothesis for where 681 came from.
All from the paper's own equations (10),(11),(13).
"""
import numpy as np

Q, EPS0 = 1.602176634e-19, 8.8541878128e-14
EPS_S = 11.7

def E_SC_Vum(tS_um):                       # eq (10), returned in V/um
    L = np.log(tS_um)
    return 9.783 * (21.765 - L) / (L + 3.975)

def E_I_Vum(tS, epsI, sig):                # eq (11)
    return (EPS_S * E_SC_Vum(tS) + Q * sig / EPS0 * 1e-4) / epsI

def V_BV(tS, epsI, sig, tI=2.0):           # eq (13)
    return (E_SC_Vum(tS) * (EPS_S * tI + 0.5 * epsI * tS)
            + Q * sig * (tI * 1e-4) / EPS0) / epsI

tS, sig = 1.0, 4e12
esc = E_SC_Vum(tS)
charge_39 = Q * sig / (3.9 * EPS0) * 1e-4      # V/um
charge_20 = Q * sig / (2.0 * EPS0) * 1e-4

print("Paper sentence under test:")
print('  "t_S=1um, eps_I=3.9eps0:  E_I = 161 (sig=0)  and  681 V/um (sig=4e12)"')
print()
print(f"E_S,C(1um) from eq(10)                    = {esc:7.2f} V/um")
print()
print("CHECK 1 - direct from eq (11) with eps_I = 3.9 eps0:")
print(f"  E_I(sig=0)      = 11.7*{esc:.2f}/3.9        = {E_I_Vum(tS,3.9,0):7.1f} V/um   (paper: 161)  OK")
print(f"  charge term     = q*4e12/(3.9*eps0)      = {charge_39:7.1f} V/um")
print(f"  E_I(sig=4e12)   = 161 + {charge_39:.0f}            = {E_I_Vum(tS,3.9,sig):7.1f} V/um   (paper: 681)  <-- MISMATCH")
print()
print("CHECK 2 - back-implied from the paper's OWN V_B,V numbers (eq 13):")
v0, v4 = V_BV(tS,3.9,0), V_BV(tS,3.9,sig)
print(f"  V_B,V(sig=0)    = {v0:6.1f} V   (paper: 348)  OK")
print(f"  V_B,V(sig=4e12) = {v4:6.1f} V   (paper: 719)  OK")
print(f"  charge dV = {v4-v0:.1f} V over t_I=2um  ->  dE_I = {(v4-v0)/2:.1f} V/um")
print(f"  ->  E_I implied by their own 719 V = 161 + {(v4-v0)/2:.0f} = {161+(v4-v0)/2:.0f} V/um  (not 681)")
print()
print("CHECK 3 - the paper's own Fig 3(a): at t_S=1, the sigma=4 curve reads ~35e5 V/cm")
print(f"  = ~350 V/um, consistent with {E_I_Vum(tS,3.9,sig):.0f}, NOT 681 (=68e5, near the sig=10 curve).")
print()
print("WHERE 681 LIKELY CAME FROM:")
print(f"  Using the LOW-K value eps_I = 2.0 eps0 instead of 3.9:")
print(f"  E_I = 11.7*{esc:.2f}/2.0 + q*4e12/(2.0*eps0)")
print(f"      = {11.7*esc/2.0:6.1f} + {charge_20:6.1f} = {E_I_Vum(tS,2.0,sig):6.1f} V/um")
print(f"  and with the paper's own rounded E_I(eps=2.0,sig=0)=319:  319 + {charge_20:.0f} = {319+charge_20:.0f}  ~= 681 !")
print()
sig_needed = (681 - E_I_Vum(tS,3.9,0)) * 1e4 * 3.9 * EPS0 / Q
print(f"  (Alternative: to get 681 at eps=3.9 you would need sigma_in = {sig_needed:.2e} cm^-2, ~11e12, not 4e12.)")
