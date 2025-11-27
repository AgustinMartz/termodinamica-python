#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Programa para hacer calculos con la ecuación de estado de van der Waals para sustancia puras

R = 8.314472  # J /mol · K

print(' ECUACIÓN DE ESTADO DE  VAN  DER  WAALS ')

TC = float(input('Temperatura crítica (en K): '))
PC = float(input('Presión crítica (en bar): '))
PC = PC * 100000

a = ((27/64) * (R * TC)**2) / (PC)
a = round(a, 5)
b = (R * TC) / (8 * PC)

print(' a = ', a, ' [(J m^3) / mol^2] ')
print(' b = ', format(b, '0.5e'),' [m^3 / mol] \n ')

tipo = input('''Elija la varible que desea calcular :
 Teclear P si va a calcular la presión
 Teclear T si va a calcular la temperatura
 Teclear v si va a calcular el volumen específico \n''')

if tipo.lower() == 'p':
    T = float(input('Temperatura (en Kelvin): '))
    v = float(input('Volumen específico (en m^3/mol): '))
    P = ((R * T)/(v - b)) - (a / v**2)
    P = round(P, 2)
    print(' La presión es igual a ', P, ' Pa   =  ', round(P/100000, 3), ' bar ')

elif tipo.lower() == 't':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    v = float(input('Volumen específico (en m^3/mol): '))
    T = (P + (a/v**2)) * ((v - b)/R)
    T = round(T, 2)
    print(' La temperatura es igual a ', T, ' K ')

elif tipo.lower() == 'v':
    P = float(input('Presión (en bar): '))
    P = P * 100000
    T = float(input('Temperatura (en Kelvin): '))

    import numpy as np
    c1 = P
    c2 = -((R*T)+(P*b))
    c3 = a
    c4 = -a*b
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    v1 = print('v1 = ', np.round(v[0], 6), '[m^3/mol] ')
    v2 = print('v2 = ', np.round(v[1], 8), '[m^3/mol] ')
    v3 = print('v3 = ', np.round(v[2], 8), '[m^3/mol] ')
    
else:
    print(' No eligió una variable válida')

# Fin
