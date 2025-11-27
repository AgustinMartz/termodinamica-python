
def vanwaals():
    import numpy as np
    import chemicals as cm
    import pyromat as pm
    import math
    import matplotlib as plt
    print(' ECUACIÓ“N DE ESTADO DE  VAN  DER  WAALS ')

    R = 8.314472

    sus = input("Ingresa la fórmula de una de las siguientes sustancias: \n "
                "Agua - H2O \n"
                "Oxígeno - O2 \n "
                "Nitrógeno - N2 \n"
                "Dióxido de carbono - CO2 \n"
                "Metano - CH4 \n"
                "R-134a - C2H2F4 \n"
                "R-1234ze(E) - C2H2F4_1 \n"
                " Otra: Ingresa 1 \n"
                ":")
    if sus == "1":
        susta = input(f'Ingresa la fórmula de la sustancia')
        sustancia = cm.search_chemical(susta)
        PC = cm.critical.Pc(sustancia.CAS)
        TC = cm.critical.Tc(sustancia.CAS)
        w = cm.acentric.omega(sustancia.CAS)
    else:
        sustancia = pm.get(f'mp.{sus}')

        TC = sustancia.critical()[0]
        PC = sustancia.critical()[1]
        PC = PC * 100000
        t = 0.7 * TC
        p = sustancia.ps(T=t)
        w = -1 - math.log(p / sustancia.critical()[1], 10)
    a = ((27 / 64) * (R * TC) ** 2) / (PC)
    a = round(a, 5)
    b = (R * TC) / (8 * PC)

    print(' a = ', a, ' [(J m^3) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] \n ')

    tipo = input('''Elija la varible que desea calcular :
     Teclear P si va a calcular la presión
     Teclear T si va a calcular la temperatura
     Teclear v si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        v = float(input('Volumen específico (en m^3/mol): '))
        P = ((R * T) / (v - b)) - (a / v ** 2)
        P = round(P, 2)
        print(' La presIón es igual a ', P, ' Pa   =  ', round(P / 100000, 3), ' bar ')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): '))
        P = P * 100000
        v = float(input('Volumen específico (en m^3/mol): '))
        T = (P + (a / v ** 2)) * ((v - b) / R)
        T = round(T, 2)
        print(' La temperatura es igual a ', T, ' K ')

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): '))
        P = P * 100000
        T = float(input('Temperatura (en Kelvin): '))


        c1 = P
        c2 = -((R * T) + (P * b))
        c3 = a
        c4 = -a * b
        coeff = [c1, c2, c3, c4]
        v = np.roots(coeff)
        vs = []

        def v_real(v):
            for i in v:
                if np.isreal(i) == True:
                    vs.append(np.real(i))

                if len(vs) == 1:
                    vreal = round(vs[0], 8)

                else:
                    print('\n Se encontraron varias raíces reales, elije una: \n')
                    for i in range(len(vs)):
                        print(f'v{i + 1} = ', np.round(v[i], 8), '[m^3/mol] \n')

                        op = int(input())
                        vreal = round(v[op - 1], 8)
            return vreal


    else:
        print(' No eligió una variable válida')

    volumen = v_real(v)

    print(f'\nEl volumen específico es igual a {volumen} [m^3/mol] ')

    lnp = (b/(volumen - b)) - np.log(((volumen-b)*P)/(R*T)) - (2*a)/(R*T*volumen)
    fug = np.exp(lnp)
    print(f'Coeficiente de fugacidad : {fug}')
    print(f'Fugacidad : {round((fug * P ) / 1000, 3)} kPa ó {round((fug*P)/100000, 3)} bar')



