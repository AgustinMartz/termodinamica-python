def redlich_mezcla():
    import chemicals as cm
    from scipy.optimize import newton
    import numpy as np
    print('REGLA DE MEZCLADO DE VAN DER WAALS:')
    print('para la ecuación de estado de Redlich-Kwong\n')

    R = 8.314472  # J /molˆ™ K
    Tcs = []
    Pcs = []
    ws = []
    ass = []
    bs = []
    ys = []

    def a_b():
        susta = input(f'Ingresa la fórmula de la sustancia : ')
        sustancia = cm.search_chemical(susta)
        Pc = cm.critical.Pc(sustancia.CAS)
        Pcs.append(Pc)
        Tc = cm.critical.Tc(sustancia.CAS)
        Tcs.append(Tc)
        w = cm.acentric.omega(sustancia.CAS)
        ws.append(w)
        a = (9*((2**(1/3))-1))**(-1) * ((R**2 * Tc**2.5)/Pc)
        ass.append(a)
        b = (((2**(1/3))-1)/3 ) * ((R * Tc) / Pc)
        bs.append(b)
        return a, b
    n = int(input('Número de especies: '))

    A = []
    B = []
    Y = []
    i = 0
    comprobacion = 0
    R = 8.314472  # J /molˆ™ K
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
    for i in range(n):
        for j in range(n):
            amix = amix + (Y[i] * Y[j] * (A[i] * A[j]) ** 0.5)
    print('amix = ', round(amix, 8), ' [(J m^3 K^0.5) / mol^2] ')

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
    c1 = 1
    c2 = -(R * T) / P
    c3 = (a / (P * T ** 0.5)) - (b ** 2) - ((b * R * T) / P)
    c4 = -(a * b) / (P * T ** 0.5)
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

    volumen1 = v_real(v)

    print(f'\nEl volumen específico es igual a {volumen1} [m^3/mol] ')

    print('\n Comparando con valores pseudocríticos:')

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

    print(f'El volumen específico es igual a {volumen} [m^3/mol] ')

    print(' \n CÁLCULO DE FUGACIDADES')

    if n == 2:
        z = (P * volumen1) / (R * T)
        G = (bs[0]/bmix)*(z - 1)
        H = np.log(((volumen1 - bmix) * P) / (R * T))
        J = 1/(bmix*R*(T**1.5))
        K = ((amix*bs[0]/bmix) - 2*(Y[1]*ass[0] + Y[0]*((ass[0]*ass[1])**0.5)))
        L = np.log( 1 + (bs[0]/volumen1))
        lnp = G - H + (J*K*L)
        fug = np.exp(lnp)

        print(f'Coeficiente de fugacidad especie 1: {fug}')
        print(f'Fugacidad de la especie 1 : {round((fug * P * Y[1]) / 100000, 3)} kPa ó {round((fug * P * Y[1]) / 100000, 3)} bar')

        z = (P * volumen1) / (R * T)
        G = (bs[1] / bmix) * (z - 1)
        H = np.log(((volumen1 - bmix) * P) / (R * T))
        J = 1 / (bmix * R * (T ** 1.5))
        K = ((amix*bs[1]/bmix) - 2*(Y[0]*ass[1] + Y[0]*((ass[1]*ass[0])**0.5)))
        L = np.log(1 + (bs[1] / volumen1))
        lnp = G - H + (J * K * L)
        fug = np.exp(lnp)

        print(f'Coeficiente de fugacidad especie 2: {fug}')
        print(f'Fugacidad de la especie 2 : {round((fug * P * Y[0]) / 1000, 3)} kPa ó {round((fug * P * Y[0]) / 100000, 3)} bar')

    elif n > 2:
        CFs = []
        Fugs = []
        Am = ass.copy()
        Ym = ys.copy()
        for j in range(n):
            mix = 0
            for i in range(n):
                mix = mix + ys[i] * (ass[j] * Am[i])**0.5

            z = (P * volumen1) / (R * T)
            G = (bs[j] / bmix) * (z - 1)
            H = np.log(((volumen1 - bmix) * P) / (R * T))
            J = 1 / (bmix * R * (T ** 1.5))
            K = ((amix*bs[j])/(bmix)) - 2*mix
            L = np.log(1 + (bmix/volumen1))
            lnp = G - H + (J*K*L)
            fug = np.exp(lnp)
            CFs.append(fug)
            Fugs.append(fug * P * ys[j])

        Am.append(Am[0])
        del Am[0]
        Ym.append(Ym[0])
        del Ym[0]

        for j in range(n):
            print(f'Coeficiente de fugacidad de la especie {j + 1} : {CFs[j]}')
            print(f'Fugacidad de la especie {j + 1} : {round(Fugs[j]/1000, 3)} kPa ó {round((Fugs[j]) / 100000, 3)} bar')

