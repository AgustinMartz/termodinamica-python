"""
vdw_pura.py
Ecuación de Van der Waals para sustancias puras

Autores:
- Estudiante de ejemplo (Universidad de Ejemplo)
Asesor:
- Dr. Agustín Martínez Ruvalcaba
Programa Delfín 2023-2024

Descripción:
Este script resuelve el volumen molar de una sustancia pura usando la ecuación de Van der Waals.
"""

import numpy as np
from scipy.optimize import fsolve

# Constantes de Van der Waals para dióxido de carbono (CO2)
a = 3.59  # L^2·bar/mol^2
b = 0.0427  # L/mol
R = 0.08314  # L·bar/K·mol

# Condiciones
T = 300  # K
P = 10   # bar

# Ecuación de Van der Waals
def vdw_eq(V):
    return (R * T) / (V - b) - a / V**2 - P

# Estimación inicial
V0 = R * T / P

# Solución
V = fsolve(vdw_eq, V0)[0]

print(f"Volumen molar estimado a T = {T} K y P = {P} bar: {V:.4f} L/mol")
