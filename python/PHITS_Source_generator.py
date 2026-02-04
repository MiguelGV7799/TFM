# PHITSFusion © 2026 by Miguel García, Manuel Cotelo
# Licensed under CC BY-NC-SA 4.0
# To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/

import argparse
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------
# Constantes / parámetros
# -------------------------
cphys = argparse.Namespace(
    Na  = 6.02214076e23,   # Avogadro [1/mol]
    awD = 2.01410178,      # D  [g/mol]
    awT = 3.01604928,      # T  [g/mol]
)

params = {
    "xD": 0.5,
    "xT": 0.5,
    "particle": "neutron",
    "kinetic_energy": 14.1,  # MeV
}

# -------------------------
# Hotspot Patel (perfil T(r), rho(r))
# -------------------------
class HotSpotPatel:

    def __init__(
        self,
        T0_keV  = 6,
        P_Gbar  = 337.0,
        R0_cm   = 3.47e-3,   # radio geométrico (PHITS)
        Rhs_cm  = 3.18e-3,   # radio de cálculo (98% yield), ajustar al shot
        b_shape = 0.67
    ):
        self.T0_keV  = float(T0_keV)
        self.P_Gbar  = float(P_Gbar)
        self.R0_cm   = float(R0_cm)
        self.Rhs_cm  = float(Rhs_cm)
        self.b_shape = float(b_shape)

        if not (0.0 < self.Rhs_cm <= self.R0_cm):
            raise ValueError("Debe cumplirse 0 < Rhs_cm <= R0_cm")

        # Densidad central isobárica (Patel): rho[g/cc] = 1.3 * P[Gbar] / T[keV]
        self.rho0_g_cc = 1.3 * self.P_Gbar / self.T0_keV

        print("Hotspot inicializado:")
        print(f"  R0  = {self.R0_cm*1e4:.2f} µm (geometría PHITS)")
        print(f"  Rhs = {self.Rhs_cm*1e4:.2f} µm (cálculos)")
        print(f"  rho0 = {self.rho0_g_cc:.2f} g/cm³")

    def __call__(self, radius_in_cm):
        r = np.asarray(radius_in_cm, dtype=float)

        # Dominio físico del perfil: 0 <= r < R0
        r_eff = np.minimum(r, self.R0_cm * (1.0 - 1e-12))
        x = r_eff / self.R0_cm  # 0 <= x < 1

        expo = 1.0 / (1.0 + self.b_shape)
        tem_in_keV = self.T0_keV * (1.0 - x**2)**expo

        # Guardarraíl numérico 
        T_min_keV = 0.1
        tem_in_keV = np.maximum(tem_in_keV, T_min_keV)

        den_in_g_cm3 = 1.3 * self.P_Gbar / tem_in_keV
        return den_in_g_cm3, tem_in_keV

    def rho_avg_up_to_Rhs(self, npts=2000):
        rs = np.linspace(0.0, self.Rhs_cm, npts)
        dens, _ = self(rs)
        integral = np.trapezoid(dens * rs**2, rs)
        return 3.0 * integral / (self.Rhs_cm**3)

# Instancia
hs = HotSpotPatel(
    T0_keV=6,
    P_Gbar=337.0,
    R0_cm=3.47e-3,
    Rhs_cm=3.18e-3,
    b_shape=0.67
)

# -------------------------
# Diagnóstico perfiles (hasta R0 para ver forma; promedios hasta Rhs)
# -------------------------
rs_plot = np.linspace(0.0, hs.R0_cm, 2000)
dens_plot, tems_plot = hs(rs_plot)

r_norm = rs_plot / hs.R0_cm
T_norm = tems_plot / hs.T0_keV
rho_norm = dens_plot / hs.rho0_g_cc

rho_avg = hs.rho_avg_up_to_Rhs()
print(f"Densidad media del hotspot hasta Rhs (g/cm^3): {rho_avg:.2f}")

# Plots
fig1, ax1 = plt.subplots(figsize=(6, 4))
ax1.plot(r_norm, T_norm, label="T/T₀")
ax1.axvline(hs.Rhs_cm/hs.R0_cm, linestyle="--", label="Rhs/R0")
ax1.set_xlabel("r / R₀")
ax1.set_ylabel("T / T₀")
ax1.set_title("Perfil de temperatura normalizado (Patel)")
ax1.grid(True, which="both", linestyle="--", alpha=0.5)
ax1.legend()
plt.show()

fig2, ax2 = plt.subplots(figsize=(6, 4))
ax2.plot(r_norm, rho_norm, label="ρ/ρ₀")
ax2.axvline(hs.Rhs_cm/hs.R0_cm, linestyle="--", label="Rhs/R0")
ax2.set_xlabel("r / R₀")
ax2.set_ylabel("ρ / ρ₀")
ax2.set_title("Perfil de densidad normalizado (Patel)")
ax2.grid(True, which="both", linestyle="--", alpha=0.5)
ax2.legend()
plt.show()

# -------------------------
# Reactividad Bosch-Hale (D+T)
# -------------------------
class Reactivity:
    def __init__(self):
        # Bosch–Hale D(T,n)4He
        self.C1 = 1.17302e-09
        self.C2 = 1.51361e-02
        self.C3 = 7.51886e-02
        self.C4 = 4.60643e-03
        self.C5 = 1.35000e-02
        self.C6 = -1.06750e-04
        self.C7 = 1.36600e-05

        self.mrc2_keV = 1124656
        self.BG = 34.3827

    def __call__(self, T_keV):
        T = np.array(T_keV, dtype=float)
        T = np.maximum(T, 1e-6)
        T = np.clip(T, 0.2, 100.0)

        C2, C3 = self.C2, self.C3
        C4, C5 = self.C4, self.C5
        C6, C7 = self.C6, self.C7

        num = T * (C2 + T * (C4 + T * C6))
        den = 1 + T * (C3 + T * (C5 + T * C7))
        theta = T / (1 - (num / den))

        xi = (self.BG**2 / (4.0 * theta))**(1/3)

        prefactor = self.C1 * theta * np.exp(-3 * xi)
        sqrt_term = np.sqrt(xi / (self.mrc2_keV * T**3))
        sv = prefactor * sqrt_term
        return sv

sv = Reactivity()
for T in [1, 5, 10, 20, 30, 40]:  #Prueba
    print(f"T = {T:4.1f} keV -> <σv> = {sv(T):.3e} cm^3/s")
# -------------------------
# Tasa local de reacciones
# -------------------------
class Rates:
    def __init__(self, func_state, func_reactivity):
        Abar = params["xD"]*cphys.awD + params["xT"]*cphys.awT  # g/mol (masa molar media)
        self.den_to_ni = cphys.Na / Abar  # [1/(g)] -> ni = rho[g/cc]*den_to_ni => [1/cc]
        self.func_state = func_state
        self.func_reactivity = func_reactivity

    def __call__(self, radius_in_cm):
        den_gcc, tem_keV = self.func_state(radius_in_cm)
        ni_cm3 = self.den_to_ni * den_gcc
        return self.func_reactivity(tem_keV) * 0.25 * ni_cm3**2

rates = Rates(hs, sv)

# -------------------------
# Discretización e integración hasta Rhs (cálculos)
# -------------------------
n_shells = 10
radius_edges = np.linspace(0.0, hs.Rhs_cm, n_shells + 1)  # [cm]
radius_centers = 0.5 * (radius_edges[1:] + radius_edges[:-1])

rows = []
for r_lo, r_hi, r_center in zip(radius_edges[:-1], radius_edges[1:], radius_centers):
    volume = (4.0/3.0) * np.pi * (r_hi**3 - r_lo**3)

    # Integración volumétrica: ∫ 4π r^2 R(r) dr
    integrand = lambda r: 4.0 * np.pi * r*r * rates(r)
    total_rate, _ = sp.integrate.quad(integrand, r_lo, r_hi, limit=200)

    rate_mean = total_rate / volume

    rows.append([r_lo, r_hi, r_center, volume, total_rate, rate_mean])

data = pd.DataFrame(
    rows,
    columns=[
        "r_lo_cm",
        "r_hi_cm",
        "r_center_cm",
        "volume_cm3",
        "rate_total_per_s",
        "rate_mean_cm3_per_s"
    ]
)

# Probabilidad por capa
data["P"] = data["rate_total_per_s"] / data["rate_total_per_s"].sum()
params["profile"] = data

print(params["profile"][["r_lo_cm","r_hi_cm","rate_total_per_s","P"]])
print("Suma P:", data["P"].sum())

# Perfil radial de tasa media
plt.figure()
plt.plot(data["r_center_cm"], data["rate_mean_cm3_per_s"], marker="o")
plt.xlabel("Radio [cm]")
plt.ylabel("Tasa volumétrica media [reacciones/(cm$^3$·s)]")
plt.title("Perfil radial de tasa volumétrica de reacciones (hasta Rhs)")
plt.grid(True)
plt.show()

# -------------------------
# Generación del bloque [Source] PHITS
# -------------------------
def phits_input_common():
    return '''
# PHITS input common settings
[ Title ]
  Hot Spot Neutron Source
'''

def phits_input_source(row, particle, kinetic_energy):
    # IMPORTANTE: 'P' es la probabilidad de la capa (row["P"])
    return f'''
<source> = {row["P"]:.8e}
  proj = {particle}
  s-type = 9 # spherical shell
  r1 = {row["r_lo_cm"]:13.5e} # [cm]
  r2 = {row["r_hi_cm"]:13.5e} # [cm]
  dir = all # isotropic
  e0 = {kinetic_energy} # [MeV]
'''

def phits_input_sources(profile, particle, kinetic_energy):
    s = "[ Source ]\n"
    for _, row in profile.iterrows():
        s += phits_input_source(row, particle, kinetic_energy)
    return s

def phits_input(params):
    return phits_input_common() + phits_input_sources(
        params["profile"],
        params["particle"],
        params["kinetic_energy"]
    )

input_text = phits_input(params)
print(input_text)

