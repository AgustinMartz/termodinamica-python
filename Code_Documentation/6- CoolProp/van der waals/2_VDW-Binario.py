import CoolProp.CoolProp as CP
from scipy.optimize import newton
import numpy as np

# Constante universal de los gases en J/(mol*K)
R = 8.3144621

print('ECUACIÓN DE ESTADO DE  VAN  DER  WAALS - Mezcla Binaria ')

# Seleccionar las sustancias en inglés
sustancia1 = input('Escribe el nombre de la primer sustancia (en inglés): ')
sustancia2 = input('Escribe el nombre de la segunda sustancia (en inglés): ')

# Obtener las fracciones molares
x1 = float(input('Escribe la fracción molar de la 1era sustancia: '))
x2 = 1 - x1

# Obtener las propiedades críticas de las sustancias
try:
    TC1 = CP.PropsSI('Tcrit', sustancia1)
    PC1 = CP.PropsSI('Pcrit', sustancia1)
    TC2 = CP.PropsSI('Tcrit', sustancia2)
    PC2 = CP.PropsSI('Pcrit', sustancia2)

    print('\n')
    print(f'1era Sustancia seleccionada: {sustancia1}')
    print(f'Temperatura crítica (TC): {TC1} K')
    print(f'Presión crítica (PC): {PC1} Pa \n')

    print(f'2da Sustancia seleccionada: {sustancia2}')
    print(f'Temperatura crítica (TC): {TC2} K')
    print(f'Presión crítica (PC): {PC2} Pa \n')

# Calcular las constantes de Van der Waals
    a1 = ((27 / 64) * (R * TC1)**2) / PC1
    b1 = (R * TC1) / (8 * PC1)
    a2 = ((27 / 64) * (R * TC2) ** 2) / PC2
    b2 = (R * TC2) / (8 * PC2)

# Constantes de mezcla
    a12 = np.sqrt(a1*a2)
    amix = x1**2 * a1 + 2 * x1 * x2 * a12 + x2**2 * a2
    bmix = x1 * b1 + x2 * b2

    print('amix = ', amix, '[(J m^3 / mol^2]')
    print('bmix = ', format(bmix, '0.5e'), '[m^3 / mol] \n ')

# Seleccionar la variable a calcular
    tipo = input('''Elija la variable que desea calcular :
    Teclear P si va a calcular la presión
    Teclear T si va a calcular la temperatura
    Teclear V si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        V = float(input('Volumen específico (en m^3/mol): '))
        P = ((R * T) / (V - bmix)) - (amix / V**2)
        P = round(P, 2)
        print(' La presión es igual a ', P, ' Pa   =  ', round(P / 100000, 3), ' bar ')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): '))
        P = P * 100000  # Convertir bar a Pa
        V = float(input('Volumen específico (en m^3/mol): '))
        T = (P + (amix / V**2)) * ((V - bmix) / R)
        T = round(T, 2)
        print(' La temperatura es igual a ', T, ' K ')

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): '))
        P = P * 100000  # Convertir bar a Pa
        T = float(input('Temperatura (en Kelvin): '))

        # Resolver la ecuación cúbica para V usando numpy.roots
        coeff = [P, -(R * T + P * bmix), amix, -amix * bmix]
        V_roots = np.roots(coeff)
        V_roots = np.real(V_roots[np.isreal(V_roots)])  # Filtrar raíces reales

        print('Volúmenes específicos posibles:')
        for i, v in enumerate(V_roots, start=1):
            print(f'v{i} = {np.round(v, 6)} [m^3/mol]')

    else:
        print(' No eligió una variable válida')
except Exception as e:
    print(f'Error al obtener las propiedades de la sustancia: {e}')

