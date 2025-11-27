#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Programa para hacer calculos con la ecuación de estado de Redlich-Kwong para sustancia puras

R = 8.314472  # J /mol · K

print(' ECUACIÓN DE ESTADO DE  REDLICH-KWONG ')

TC = float(input('Temperatura crítica (en K): '))
PC = float(input('Presión crítica (en bar): '))
PC = PC * 100000

a = (9*((2**(1/3))-1))**(-1) * ((R**2 * TC**2.5)/PC)
a = round(a, 5)
b = (((2**(1/3))-1)/3 ) * ((R * TC) / PC)
print(' a = ', a, ' [(J m^3 K^0.5) / mol^2] ')
print(' b = ', format(b, '0.5e'),' [m^3 / mol] \n ')

tipo = input('''Elija la varible que desea calcular :
 Teclear P si va a calcular la presión
 Teclear T si va a calcular la temperatura
 Teclear v si va a calcular el volumen específico \n''')

if tipo.lower() == 'p':
    T = float(input('Temperatura (en Kelvin): '))
    v = float(input('Volumen específico (en m^3/mol): '))
    P = ((R * T)/(v - b)) - (a / (T**0.5 * v * (v + b)))
    P = round(P, 2)
    print(' La presión es igual a ', P, ' Pa   =  ', round(P/100000, 3), ' bar ')

elif tipo.lower() == 't':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    v = float(input('Volumen específico (en m^3/mol): '))
    m1 = -(P * (v - b)) / R
    m2 = -a / (R * v)

    from scipy.optimize import newton
    def f(T):
        return T**(3/2) + (m1 * T**0.5) + m2
    T = TC
    T0 = newton(f, T, fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
    T0 = round(T0.real, 2)
    print(' La temperatura es igual a ', T0, ' K ')

elif tipo.lower() == 'v':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    T = float(input('Temperatura (en Kelvin): '))

    import numpy as np
    c1 = 1
    c2 = -(R*T)/P
    c3 = (a / (P * T**0.5)) - (b**2) - ((b * R * T) / P)
    c4 = -(a*b)/(P * T**0.5)
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    v1 = print('v1 = ', np.round(v[0], 6), '[m^3/mol] ')
    v2 = print('v2 = ', np.round(v[1], 8), '[m^3/mol] ')
    v3 = print('v3 = ', np.round(v[2], 8), '[m^3/mol] ')
    
else:
    print(' No eligió una variable válida')

# Fin

