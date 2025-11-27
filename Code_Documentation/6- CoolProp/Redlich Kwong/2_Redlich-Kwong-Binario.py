# ECUACIÓN DE ESTADO DE REDLICH--KWONG PARA UNA MEZCLA BINARIA
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import newton

# Constante de los gases
R = 8.314472  # J /mol K

print('-----ECUACIÓN DE ESTADO DE REDLICH-KWONG-----')

# Selecciona la sustancia en la base de CoolProb (en ingles)
nombre_sustancia1 = input('Ingrese el nombre de la primera sustancia en ingles: ')
nombre_sustancia2 = input('Ingrese el nombre de la segunda sustancia en ingles: ')
print("")

# Obtener propiedades críticas de cada sustancia usando CoolProp
try:
    TC1 = CP.PropsSI('Tcrit', nombre_sustancia1)
    PC1 = CP.PropsSI('Pcrit', nombre_sustancia1) * 1e-5  # Convertir Pa a bar
    print(f'Sustancia seleccionada: {nombre_sustancia1}')
    print(f'Temperatura crítica (TC): {TC1} K')
    print(f'Presión crítica (PC): {PC1} bar')
    print("")
    TC2 = CP.PropsSI('Tcrit', nombre_sustancia2)
    PC2 = CP.PropsSI('Pcrit', nombre_sustancia2) * 1e-5  # Convertir Pa a bar
    print(f'Sustancia seleccionada: {nombre_sustancia2}')
    print(f'Temperatura crítica (TC): {TC2} K')
    print(f'Presión crítica (PC): {PC2} bar')
except Exception as e:
    print(f'Error: {e}')
    print('Sustancia no encontrada en la base de datos de CoolProp.')
    exit()

# Fracciones molares
print("")
x1 = float(input('Ingrese la fracción molar de la primera sustancia: '))
x2 = 1 - x1

# Parámetros de cada sustancia de la ecuación de Redlich-Kwong
a1 = (9*((2**(1/3))-1))**(-1) * ((R**2 * TC1**2.5)/(PC1*1e5))
a1 = round(a1, 5)
b1 = (((2**(1/3))-1)/3) * ((R * TC1) / (PC1*1e5))
a2 = (9*((2**(1/3))-1))**(-1) * ((R**2 * TC2**2.5)/(PC2*1e5))
a2 = round(a2, 5)
b2 = (((2**(1/3))-1)/3) * ((R * TC2) / (PC2*1e5))

# Parámetros de la mezcla utilizando la regla de mezclado
print("")
a = x1**2 * a1 + x2**2 * a2 + 2 * x1 * x2 * (a1 * a2)**0.5
b = x1 * b1 + x2 * b2
print(' amix = ', a, ' [(J m^3 K^0.5) / mol^2] ')
print(' bmix = ', format(b, '0.5e'), ' [m^3 / mol] \n ')

# Elección de la variable a calcular
tipo = input('''Elija la variable que desea calcular :
 Teclear P si va a calcular la presión
 Teclear T si va a calcular la temperatura
 Teclear v si va a calcular el volumen específico \n''')

# Cálculo de presión
if tipo.lower() == 'p':
    T = float(input('Temperatura (en Kelvin): '))
    v = float(input('Volumen específico (en m^3/mol): '))
    P = ((R * T)/(v - b)) - (a / (T**0.5 * v * (v + b)))
    P = round(P, 2)
    print(' La presión es igual a ', P, ' Pa   =  ', round(P/100000, 3), ' bar ')

# Cálculo de temperatura
elif tipo.lower() == 't':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    v = float(input('Volumen específico (en m^3/mol): '))
    m1 = -(P * (v - b)) / R
    m2 = -a / (R * v)

    def f(T):
        return T**(3/2) + (m1 * T**0.5) + m2
    T0 = newton(f, max(TC1, TC2), fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
    T0 = round(T0.real, 2)
    print(' La temperatura es igual a ', T0, ' K ')

# Cálculo de volumen
elif tipo.lower() == 'v':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    T = float(input('Temperatura (en Kelvin): '))

    c1 = 1
    c2 = -(R*T)/P
    c3 = (a / (P * T**0.5)) - (b**2) - ((b * R * T) / P)
    c4 = -(a*b)/(P * T**0.5)
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    print('v1 = ', np.round(v[0], 6), '[m^3/mol] ')
    print('v2 = ', np.round(v[1], 8), '[m^3/mol] ')
    print('v3 = ', np.round(v[2], 8), '[m^3/mol] ')

else:
    print(' No eligió una variable válida')