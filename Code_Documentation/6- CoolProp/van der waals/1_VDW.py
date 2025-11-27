import CoolProp.CoolProp as CP
import numpy as np

# Constante universal de los gases en J/(mol*K)
R = 8.3144621

print('ECUACIÓN DE ESTADO DE VAN DER WAALS')

# Seleccionar la sustancia en inglés
sustancia = input('Escribe el nombre de la sustancia en inglés: ')

# Obtener las propiedades críticas
try:
    TC = CP.PropsSI('Tcrit', sustancia)
    PC = CP.PropsSI('Pcrit', sustancia)
    molar_mass = CP.PropsSI('M', sustancia)  # Obtener la masa molar
    print(f'Sustancia seleccionada: {sustancia}')
    print(f'Temperatura crítica (TC): {TC} K')
    print(f'Presión crítica (PC): {PC} Pa')
    print(f'Masa molar: {molar_mass} g/mol')

    # Calcular las constantes de Van der Waals
    a = ((27 / 64) * (R * TC)**2) / PC
    b = (R * TC) / (8 * PC)

    print('a = ', a, ' [(J m^3) / mol^2] ')
    print('b = ', format(b, '0.5e'), ' [m^3 / mol] \n ')

    # Seleccionar la variable a calcular
    tipo = input('''Elija la variable que desea calcular :
    Teclear P si va a calcular la presión
    Teclear T si va a calcular la temperatura
    Teclear V si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        V = float(input('Volumen específico (en m^3/mol): '))
        P = ((R * T) / (V - b)) - (a / V**2)
        P = round(P, 2)
        print('La presión es igual a ', P, ' Pa   =  ', round(P / 100000, 3), ' bar')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): '))
        P = P * 100000  # Convertir bar a Pa
        V = float(input('Volumen específico (en m^3/mol): '))
        T = (P + (a / V**2)) * ((V - b) / R)
        T = round(T, 2)
        print('La temperatura es igual a ', T, ' K')

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): '))
        P = P * 100000  # Convertir bar a Pa
        T = float(input('Temperatura (en Kelvin): '))

        # Resolver la ecuación cúbica para V usando numpy.roots
        coeff = [P, -(R * T + P * b), a, -a * b]
        V_roots = np.roots(coeff)
        V_roots = np.real(V_roots[np.isreal(V_roots)])  # Filtrar raíces reales

        print('Volúmenes específicos posibles:')
        for i, v in enumerate(V_roots, start=1):
            v_m3_kg = ( v / (molar_mass * 1e3) ) * 1000  # Convertir de m^3/mol a m^3/kg
            print(f'v{i} = {np.round(v, 6)} [m^3/mol] = {np.round(v_m3_kg, 6)} [m^3/kg]')

    else:
        print('No eligió una variable válida')
except Exception as e:
    print(f'Error al obtener las propiedades de la sustancia: {e}')
