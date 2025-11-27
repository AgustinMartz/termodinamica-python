def pengMezcla():
    import chemicals as cm
    import numpy as np
    print('REGLA DE MEZCLADO DE VAN DER WAALS:')
    print('para la ecuación de estado de Peng-Robinson\n')
    R = 8.314472  # J /molˆ™ K
    Tcs = []
    Pcs = []
    ws = []
    bs = []
    ass = []
    ys = []
    alfs = []
    def a_b_alfa(T):
        susta = input(f'Ingresa la fórmula de la sustancia : ')
        sustancia = cm.search_chemical(susta)
        Pc = cm.critical.Pc(sustancia.CAS)
        Pcs.append(Pc)
        Tc = cm.critical.Tc(sustancia.CAS)
        Tcs.append(Tc)
        w = cm.acentric.omega(sustancia.CAS)
        ws.append(w)
        Tr = T / Tc
        a = 0.45724 * ((R * Tc)**2 / Pc)
        ass.append(a)
        b = 0.07780 * ((R * Tc) / Pc)
        bs.append(b)
        k = 0.37464 + (1.54226 * w) - (0.26992 * w ** 2)
        alfa = (1 + (k * (1 - Tr ** 0.5))) ** 2
        alfs.append(alfa)
        return a, b, alfa

    n = int(input('Número de especies: '))


    A = []
    B = []
    ALFA = []
    Y = []
    i = 0
    comprobacion = 0

    T1s = []

    P1s = []

    w1s = []

    T = float(input('\nTemperatura de la mezcla (en K): '))

    while i <= n - 1:
        print('\nCalculo de "a" y "b" para la especie ', i + 1, ': ')
        a, b, alfa = a_b_alfa(T)
        A.insert(0, a)
        B.insert(0, b)
        ALFA.insert(0, alfa)
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

    for i in range(n):
        for j in range(n):
            amix = amix + (Y[i] * Y[j] * (A[i] * A[j]) ** 0.5 * (ALFA[i] * ALFA[j]) ** 0.5)
    print('amix = ', round(amix, 8), ' [(J m^3) / mol^2] ')

    j = 0
    bmix = 0
    alfamix = 0
    for i in range(n):
        bmix = bmix + (Y[i] * B[i])
    print('bmix = ', format(bmix, '0.5e'), ' [m^3 / mol] ')

    for i in range(n):
        for j in range(n):
            alfamix = alfamix + Y[i] * Y[j] * ((ALFA[i] * ALFA[j]) ** 0.5)
    print(f'alfa = {alfamix}')

    a = amix
    a = round(a, 5)
    b = bmix


    P = float(input('Presión (en bar): '))
    P = P * 100000
    alfa = alfamix
    n1 = R * T
    n2 = a * alfa

    c1 = P
    c2 = (P * b) - n1
    c3 = n2 - (3 * P * b ** 2) - (2 * n1 * b)
    c4 = (P * b ** 3) + (n1 * b ** 2) - (n2 * b)
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    real_roots = np.real(v)

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

    volumen1 = v_real(v)

    print(f'El volumen específico es igual a {volumen1} [m^3/mol] ')

    print('\n Comparando con valores pseudocríticos:')

    a = 0.45724 * ((R * Tpc) ** 2 / Ppc)
    a = round(a, 5)
    b = 0.07780 * ((R * Tpc) / Ppc)
    k = 0.37464 + (1.54226 * wpc) - (0.26992 * wpc ** 2)
    k = round(k, 5)
    print(' a = ', a, ' [(J m^3) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] ')
    print(' k = ', k, '\n')

    Tr = T / Tpc
    alfa = (1 + (k * (1 - Tr ** 0.5))) ** 2
    n1 = R * T
    n2 = a * alfa

    c1 = P
    c2 = (P * b) - n1
    c3 = n2 - (3 * P * b ** 2) - (2 * n1 * b)
    c4 = (P * b ** 3) + (n1 * b ** 2) - (n2 * b)
    coeff = [c1, c2, c3, c4]
    v = np.roots(coeff)
    real_roots = np.real(v)
    volumen = None

    vs = []

    volumen = v_real(v)

    print(f' El volumen específico es igual a {volumen} [m^3/mol] ')

    print('\n CÁLCULO DE FUGACIDADES')

    if n == 2:
        z = (P * volumen1) / (R * T)
        A = (bs[0] / bmix) * (z - 1)
        B = np.log(((volumen1 - bmix) * P) / (R * T))
        C = (amix * alfamix) / (2 * 2 ** 0.5 * bmix * R * T)
        D = (bs[0] / bmix) - ((2 / (amix * alfamix)) * (Y[1] * (ass[0] * ALFA[1]) + (Y[0]) * (ass[0] * ALFA[1]*ass[1] * ALFA[0])**0.5))
        E = np.log((volumen1 + ((1 + 2 ** 0.5) * bmix)) / (volumen1 + ((1 - 2 ** 0.5) * bmix)))
        LNCF_i = A - B + (C * D * E)
        fug = np.exp(LNCF_i)
        print(f'Coeficiente de fugacidad especie 1: {fug}')
        print(f'Fugacidad de la especie 1 : {round((fug * P * Y[1]) / 100000, 3)} kPa ó {round((fug * P * Y[1]) / 100000, 3)} bar')

        z = (P * volumen1) / (R * T)
        A = (bs[1] / bmix) * (z - 1)
        B = np.log(((volumen1 - bmix) * P) / (R * T))
        C = (amix * alfamix) / (2 * 2 ** 0.5 * bmix * R * T)
        D = (bs[1] / bmix) - ((2 / (amix * alfamix)) * (Y[0] * (ass[1] * ALFA[0]) + (Y[1]) * (ass[0] * ALFA[1]*ass[1] * ALFA[0])**0.5))
        E = np.log((volumen1 + ((1 + 2 ** 0.5) * bmix)) / (volumen1 + ((1 - 2 ** 0.5) * bmix)))
        LNCF_i = A - B + (C * D * E)
        fug = np.exp(LNCF_i)

        print(f'Coeficiente de fugacidad especie 2: {fug}')
        print(f'Fugacidad de la especie 2 : {round((fug * P * Y[0]) / 100000, 3)} kPa ó {round((fug * P * Y[0]) / 100000, 3)} bar')

    elif n > 2:
        CFs = []
        Fugs = []
        Am = ass.copy()
        Ym = ys.copy()
        alfsm = alfs.copy()
        for j in range(n):
            mix = 0
            for i in range(n):
                m = (ass[j] * ass[i] * alfs[j] * alfs[i])**0.5
                mix = mix + ys[i] * (m)

            z = (P * volumen1) / (R * T)
            G = (bs[j] / bmix) * (z - 1)
            H = np.log(((volumen1 - bmix) * P) / (R * T))
            J = (amix * alfamix) / (2 * (2 ** 0.5) * bmix * R * T)
            K = (bs[j] / bmix)
            N = (2 / (amix * alfamix)) * mix
            L = np.log((volumen1 + ((1 + 2 ** 0.5) * bmix)) / (volumen1 + ((1 - 2 ** 0.5) * bmix)))
            lnp = G - H + (J*(K-N)*L)
            fug = np.exp(lnp)
            CFs.append(fug)
            Fugs.append(fug * P * ys[j])

        for j in range(n):
            print(f'Coeficiente de fugacidad de la especie {j + 1} : {CFs[j]}')
            print(f'Fugacidad de la especie {j + 1} : {round(Fugs[j] / 1000, 3)} kPa ó {round((Fugs[j]) / 100000, 3)} bar')







