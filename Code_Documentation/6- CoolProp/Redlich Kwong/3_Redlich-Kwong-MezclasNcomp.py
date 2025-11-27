# ECUACIÓN DE ESTADO DE REDLICH--KWONG PARA UNA MEZCLA MULTICOMPONENTE
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import newton

# Constante de los gases
R = 8.314472  # J /mol K

print('-----ECUACIÓN DE ESTADO DE REDLICH-KWONG-----')

# Número de componentes
n = int(input('Ingrese el número de componentes de la mezcla: '))

# Listas para almacenar nombres de sustancias y propiedades críticas
nombres_sustancias = []
TC = []
PC = []

# Obtener propiedades críticas de cada sustancia usando CoolProp
for i in range(n):
    nombre_sustancia = input(f'Ingrese el nombre de la sustancia {i+1} en ingles: ')
    nombres_sustancias.append(nombre_sustancia)
    try:
        TC_i = CP.PropsSI('Tcrit', nombre_sustancia)
        PC_i = CP.PropsSI('Pcrit', nombre_sustancia) * 1e-5  # Convertir Pa a bar
        TC.append(TC_i)
        PC.append(PC_i)
        print(f'Sustancia seleccionada: {nombre_sustancia}')
        print(f'Temperatura crítica (TC): {TC_i} K')
        print(f'Presión crítica (PC): {PC_i} bar')
        print("")
    except Exception as e:
        print(f'Error: {e}')
        print('Sustancia no encontrada en la base de datos de CoolProp.')
        exit()

# Fracciones molares
x = []
for i in range(n):
    xi = float(input(f'Ingrese la fracción molar de la sustancia {nombres_sustancias[i]}: '))
    x.append(xi)

# Verificar que la suma de las fracciones molares sea 1
if not np.isclose(sum(x), 1.0):
    print('Las fracciones molares no suman 1.')
    exit()

# Parámetros de cada sustancia de la ecuación de Redlich-Kwong
a = []
b = []
for i in range(n):
    a_i = (9*((2**(1/3))-1))**(-1) * ((R**2 * TC[i]**2.5)/(PC[i]*1e5))
    a_i = round(a_i, 5)
    b_i = (((2**(1/3))-1)/3) * ((R * TC[i]) / (PC[i]*1e5))
    a.append(a_i)
    b.append(b_i)

# Parámetros de la mezcla utilizando la regla de mezclado
print("")
a_mix = sum(x[i] * x[j] * (a[i] * a[j])**0.5 for i in range(n) for j in range(n))
b_mix = sum(x[i] * b[i] for i in range(n))
print(' amix = ', a_mix, ' [(J m^3 K^0.5) / mol^2] ')
print(' bmix = ', format(b_mix, '0.5e'), ' [m^3 / mol] ')
print("")

# Elección de la variable a calcular
tipo = input('''Elija la variable que desea calcular :
 Teclear P si va a calcular la presión
 Teclear T si va a calcular la temperatura
 Teclear v si va a calcular el volumen específico \n''')

# Cálculo de presión
if tipo.lower() == 'p':
    T = float(input('Temperatura (en Kelvin): '))
    v = float(input('Volumen específico (en m^3/mol): '))
    P = ((R * T)/(v - b_mix)) - (a_mix / (T**0.5 * v * (v + b_mix)))
    P = round(P, 2)
    print(' La presión es igual a ', P, ' Pa   =  ', round(P/100000, 3), ' bar ')

# Cálculo de temperatura
elif tipo.lower() == 't':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    v = float(input('Volumen específico (en m^3/mol): '))
    m1 = -(P * (v - b_mix)) / R
    m2 = -a_mix / (R * v)

    def f(T):
        return T**(3/2) + (m1 * T**0.5) + m2
    T0 = newton(f, max(TC), fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
    T0 = round(T0.real, 2)
    print(' La temperatura es igual a ', T0, ' K ')

# Cálculo de volumen
elif tipo.lower() == 'v':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    T = float(input('Temperatura (en Kelvin): '))

    c1 = 1
    c2 = -(R*T)/P
    c3 = (a_mix / (P * T**0.5)) - (b_mix**2) - ((b_mix * R * T) / P)
    c4 = -(a_mix*b_mix)/(P * T**0.5)
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    print('v1 = ', np.round(v[0], 6), '[m^3/mol] ')
    print('v2 = ', np.round(v[1], 8), '[m^3/mol] ')
    print('v3 = ', np.round(v[2], 8), '[m^3/mol] ')

else:
    print(' No eligió una variable válida')
