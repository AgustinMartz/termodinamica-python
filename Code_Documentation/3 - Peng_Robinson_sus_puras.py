#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Programa para hacer cálculos con la ecuación de estado de Peng-Robinson para sustancias puras

R = 8.314472  # J /mol ∙ K

print(' ECUACIÓN DE ESTADO DE  PENG-ROBINSON ')

TC = float(input('Temperatura crítica (en K): '))
PC = float(input('Presión crítica (en bar): '))
PC = PC * 100000
w = float(input('Factor acéntrico: '))

a = 0.45724 * ((R * TC)**2 / PC)
a = round(a, 5)
b = 0.07780 * ((R * TC) / PC)
k = 0.37464 + (1.54226 * w) - (0.26992 * w**2)
k = round(k, 5)
print(' a = ', a, ' [(J m^3) / mol^2] ')
print(' b = ', format(b, '0.5e'),' [m^3 / mol] ')
print(' k = ', k, '\n')

tipo = input('''Elija la varible que desea calcular :
 Teclear P si va a calcular la presión
 Teclear T si va a calcular la temperatura
 Teclear v si va a calcular el volumen específico \n''')

if tipo.lower() == 'p':
    T = float(input('Temperatura (en Kelvin): '))
    v = float(input('Volumen específico (en m^3/mol): '))
    Tr = T / TC
    alfa = (1 + (k * (1 - Tr**0.5)))**2
    P = ((R * T)/(v - b)) - ((a * alfa) / (( v * (v + b)) + (b * (v - b))))
    P = round(P, 2)
    print(' La presión es igual a ', P, ' Pa   =  ', round(P/100000, 3), ' bar ')

elif tipo.lower() == 't':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    v = float(input('Volumen específico (en m^3/mol): '))
    g1 = R / (v - b)
    g2 = a / ((v * (v + b)) + (b * (v - b)))

    from scipy.optimize import newton
    def f(T):
        return P - (g1 * T) + (g2 * (1 + (k *(1 - (T / TC)**0.5)))**2)
    T = TC
    T0 = newton(f, T, fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
    T0 = round(T0.real, 3)
    print(' La temperatura es igual a ', T0, ' K ')

elif tipo.lower() == 'v':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    T = float(input('Temperatura (en Kelvin): '))
    Tr = T / TC
    alfa = (1 + (k * (1 - Tr**0.5)))**2
    n1 = R * T
    n2 = a * alfa

    import numpy as np
    c1 = P
    c2 = (P * b) - n1
    c3 = n2 - (3 * P * b**2) - (2 * n1 * b)
    c4 = (P * b**3) + (n1 * b**2) - (n2 * b)
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    v1 = print('v1 = ', np.round(v[0], 6), '[m^3/mol] ')
    v2 = print('v2 = ', np.round(v[1], 8), '[m^3/mol] ')
    v3 = print('v3 = ', np.round(v[2], 8), '[m^3/mol] ')

else:
    print(' No eligió una variable válida')

# Fin
