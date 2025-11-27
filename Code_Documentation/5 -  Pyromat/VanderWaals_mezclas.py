import matplotlib.pyplot as plt


def vandermezcla():
    import chemicals as cm
    import numpy as np

    print('REGLA DE MEZCLADO DE VAN DER WAALS:')

    R = 8.314472  # J /molˆ™ K
    Tcs = []
    Pcs = []
    ws = []
    ass = []
    bs = []
    ys = []

    def a_b():
        susta = input(f'Ingresa la fórmula de la sustancia :')
        sustancia = cm.search_chemical(susta)
        Pc = cm.critical.Pc(sustancia.CAS)
        Pcs.append(Pc)
        Tc = cm.critical.Tc(sustancia.CAS)
        Tcs.append(Tc)
        w = cm.acentric.omega(sustancia.CAS)
        ws.append(w)
        a = ((27 / 64) * (R * Tc) ** 2) / (Pc)
        ass.append(a)
        b = (R * Tc) / (8 * Pc)
        bs.append(b)
        return a, b



    n = int(input('Número de especies: '))

    A = []
    B = []
    Y = []
    i = 0
    comprobacion = 0
    T1s = []

    P1s = []

    w1s = []

    while i <= n - 1:
        print('Calculo de "a" y "b" para la especie ', i + 1, ': ')
        a, b = a_b()
        A.insert(0, a)
        B.insert(0, b)
        print('Fracción mol de la especie ', i + 1, ': ')
        y = float(input())
        ys.append(y)
        T1s.append(y * Tcs[i])
        P1s.append(y * Pcs[i])
        w1s.append(y * ws[i])
        comprobacion = comprobacion + y
        Y.insert(0, y)
        i = i + 1
    if comprobacion == 1:
        print(' ')
    else:
        print('\nLa suma de las fracciones mol es diferente de 1')

    i = 0
    amix = 0

    Tpc = sum(T1s)
    Ppc = sum(P1s)
    wpc = sum(w1s)
    print(Y)
    for i in range(n):
        for j in range(n):
            amix = amix + (Y[i] * Y[j] * (A[i] * A[j]) ** 0.5)
    print('amix = ', round(amix, 8), ' [(J m^3) / mol^2] ')

    j = 0
    bmix = 0

    for i in range(n):
        bmix = bmix + (Y[i] * B[i])
    print('bmix = ', format(bmix, '0.5e'), ' [m^3 / mol] \n ')

    a = amix
    a = round(a, 5)
    b = bmix

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
            vreal = round(vs[op - 1], 8)
        return vreal

    volumen1 = v_real(v)

    print(f'\nEl volumen específico es igual a {volumen1} [m^3/mol] ')

    print('Comparando con valores pseudocríticos:')

    a = ((27 / 64) * (R * Tpc) ** 2) / (Ppc)
    a = round(a, 5)
    b = (R * Tpc) / (8 * Ppc)
    print(' a = ', a, ' [(J m^3) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] \n ')

    c1 = P
    c2 = -((R * T) + (P * b))
    c3 = a
    c4 = -a * b
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    vs = []

    volumen = v_real(v)

    print(f'\nEl volumen específico es igual a {volumen} [m^3/mol] ')

    print(' \n CÁLCULO DE FUGACIDADES:')


    if n == 2:
        aa = (ass[0]*ass[1])**0.5
        lnp = (bs[0] / (volumen1 - bmix)) - np.log(((volumen1 - bmix) * P) / (R * T)) - (2 * (Y[1] * ass[0] + Y[0] * aa)) / (
                    R * T * volumen1)
        fug = np.exp(lnp)
        print(f'Coeficiente de fugacidad especie 1: {fug}')
        print(f'Fugacidad de la especie 1 : {round((fug * P * Y[0]) / 1000, 3)} kPa ó {round((fug * P * Y[1]) / 100000, 3)} bar')


        lnp = (bs[1] / (volumen1 - bmix)) - np.log(((volumen1 - bmix) * P) / (R * T)) - (
                    2 * (Y[0] * ass[1] + Y[1] * aa)) / (R * T * volumen1)
        fug = np.exp(lnp)
        print(f'Coeficiente de fugacidad especie 2: {fug}')
        print(f'Fugacidad de la especie 2 : {round((fug * P * Y[1]) / 1000, 3)} kPa ó {round((fug * P * Y[0]) / 100000, 3)} bar')


    elif n > 2:
        CFs = []
        Fugs = []
        Am = ass.copy()
        Ym = ys.copy()
        for j in range(n):
            mix = 0
            for i in range(n):
                mix = mix + Ym[i] * (ass[j] * Am[i])**0.5

            lnp = (bs[j] / (volumen1 - bmix)) - np.log(((volumen1 - b) * P) / (R * T)) - ((2 * mix) / (R * T * volumen1))
            fug = np.exp(lnp)
            CFs.append(fug)
            Fugs.append(fug * P )

        Am.append(Am[0])
        del Am[0]
        Ym.append(Ym[0])
        del Ym[0]
        for j in range(n):
            print(f'Coeficiente de fugacidad de la especie {j + 1} : {CFs[j]}')
            print(f'Fugacidad de la especie {j + 1} : {round((Fugs[j]*ys[j])/1000, 3)} kPa ó {round((Fugs[j]*ys[j]) / 100000, 3)} bar')























