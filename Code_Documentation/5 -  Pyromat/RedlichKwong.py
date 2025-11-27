def redlich():
    import numpy as np
    from scipy.optimize import newton
    import pyromat as pm
    import chemicals as cm
    import math
    print(' ECUACIÓ“N DE ESTADO DE  REDLICH-KWONG ')

    R = 8.314472  # J /mol âˆ™ K

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
        susta = input(f'Ingresa la fórmula de la sustancia:')
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

    a = (9 * ((2 ** (1 / 3)) - 1)) ** (-1) * ((R ** 2 * TC ** 2.5) / PC)
    a = round(a, 5)
    b = (((2 ** (1 / 3)) - 1) / 3) * ((R * TC) / PC)
    print(' a = ', a, ' [(J m^3 K^0.5) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] \n ')

    tipo = input('''Elija la varible que desea calcular :
     Teclear P si va a calcular la presión
     Teclear T si va a calcular la temperatura
     Teclear v si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        v = float(input('Volumen específico (en m^3/mol): '))
        P = ((R * T) / (v - b)) - (a / (T ** 0.5 * v * (v + b)))
        P = round(P, 2)
        print(' La presión es igual a ', P, ' Pa   =  ', round(P / 100000, 3), ' bar ')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): '))
        P = P * 100000
        v = float(input('Volumen específico (en m^3/mol): '))
        m1 = -(P * (v - b)) / R
        m2 = -a / (R * v)


        def f(T):
            return T ** (3 / 2) + (m1 * T ** 0.5) + m2

        T = TC
        T0 = newton(f, T, fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
        T0 = round(T0.real, 2)
        print(' La temperatura es igual a ', T0, ' K ')

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): '))
        P = P * 100000
        T = float(input('Temperatura (en Kelvin): '))
        c1 = 1
        c2 = -(R * T) / P
        c3 = (a / (P * T ** 0.5)) - (b ** 2) - ((b * R * T) / P)
        c4 = -(a * b) / (P * T ** 0.5)
        coeff = [c1, c2, c3, c4]
        v = np.roots(coeff)
        v1 = print('v1 = ', np.round(v[0], 6), '[m^3/mol] ')
        v2 = print('v2 = ', np.round(v[1], 8), '[m^3/mol] ')
        v3 = print('v3 = ', np.round(v[2], 8), '[m^3/mol] ')
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
                        print(f'v{i + 1} = ', np.round(vs[i], 8), '[m^3/mol] \n')

                        op = int(input())
                        vreal = round(vs[op - 1], 8)
            return vreal
    else:
        print(' No eligió una variable válida')

    volumen = v_real(v)

    print(f'\nEl volumen específico es igual a {volumen} [m^3/mol] ')

    z = (P*volumen)/(R*T)
    G = np.log(((volumen-b)*P)/(R*T))
    H = (a/(b*R*(T**1.5)))
    J = np.log(1 + (b/volumen))
    lnp = z - 1 - G - (H*J)
    fug = np.exp(lnp)
    print(f'Coeficiente de fugacidad : {fug}')
    print(f'Fugacidad : {round((fug * P) / 1000, 3)} kPa ó {round((fug*P)/100000, 3)} bar')
    return

