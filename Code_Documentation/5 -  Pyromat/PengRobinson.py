def pengro():
    from scipy.optimize import newton
    import numpy as np
    import pyromat as pm
    import math
    import chemicals as cm

    print(' ECUACIÓN DE ESTADO DE  PENG-ROBINSON ')
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

    a = 0.45724 * ((R * TC) ** 2 / PC)
    a = round(a, 5)
    b = 0.07780 * ((R * TC) / PC)
    k = 0.37464 + (1.54226 * w) - (0.26992 * w ** 2)
    k = round(k, 5)
    print(' a = ', a, ' [(J m^3) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] ')
    print(' k = ', k, '\n')

    tipo = input('''Elija la varible que desea calcular :
     Teclear P si va a calcular la presión
     Teclear T si va a calcular la temperatura
     Teclear v si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        v = float(input('Volumen específico (en m^3/mol): '))
        Tr = T / TC
        alfa = (1 + (k * (1 - Tr ** 0.5))) ** 2
        P = ((R * T) / (v - b)) - ((a * alfa) / ((v * (v + b)) + (b * (v - b))))
        P = round(P, 2)
        print(' La presión es igual a ', P, ' Pa   =  ', round(P / 100000, 3), ' bar ')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): '))
        P = P * 100000
        v = float(input('Volumen específico (en m^3/mol): '))
        g1 = R / (v - b)
        g2 = a / ((v * (v + b)) + (b * (v - b)))

        def f(T):
            return P - (g1 * T) + (g2 * (1 + (k * (1 - (T / TC) ** 0.5))) ** 2)

        T = TC
        T0 = newton(f, T, fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
        T0 = round(T0.real, 3)
        print(' La temperatura es igual a ', T0, ' K ')
        alfa = ((R*T)/(v-b) - P)/ (a/((v**2) + (2*b*v) - b**2))

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): '))
        P = P * 100000
        T = float(input('Temperatura (en Kelvin): '))
        Tr = T / TC
        alfa = (1 + (k * (1 - Tr ** 0.5))) ** 2
        n1 = R * T
        n2 = a * alfa

        c1 = P
        c2 = (P * b) - n1
        c3 = n2 - (3 * P * b ** 2) - (2 * n1 * b)
        c4 = (P * b ** 3) + (n1 * b ** 2) - (n2 * b)
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

        volumen = v_real(v)

    print(f'\nEl volumen específico es igual a {volumen} [m^3/mol] ')


    z = (P * volumen) / (R * T)
    lnp = z - 1  - np.log(((volumen-b)*P)/(R*T)) - ((a*alfa)/(2*(2**0.5)*b*R*T))* np.log((volumen + (1 + (2**0.5))*b)/(volumen + (1 - (2**0.5))*b))
    fug = np.exp(lnp)
    print(f'Coeficiente de fugacidad : {fug}')
    print(f'Fugacidad : {round((fug * P) / 1000, 3)} kPa ó {round((fug*P)/100000,3)} bar')


