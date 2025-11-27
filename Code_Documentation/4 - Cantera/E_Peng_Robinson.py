# Principal
def peng_rob(pc, tc, w, x, volumen, esp):

    print('------ECUACIÓN DE ESTADO DE  PENG-ROBINSON------')

    if volumen == 1:
        tipo = 'v'
    else:
        if len(x) == 1:
            tipo = input('''Elija la varible que desea calcular :
            Teclear P si va a calcular la presión
            Teclear T si va a calcular la temperatura
            Teclear v si va a calcular el volumen específico \n''')
        else:
            tipo = input('''Elija la varible que desea calcular :
            Teclear P si va a calcular la presión
            Teclear v si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        v = float(input('Volumen específico (en m^3/mol): '))
        press(tc, pc, w, x, esp, T, v)

        if len(x) == 2:
            fug_bin(tc, pc, w, esp, T)

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): ')) * 1e5
        v = float(input('Volumen específico (en m^3/mol): '))
        temp(tc, pc, w, P, v)

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): ')) * 1e5
        T = float(input('Temperatura (en Kelvin): '))
        vol(tc, pc, w, x, esp, P, T)

        if len(x) == 2:
            fug_bin(tc, pc, w, esp, T)

    else:
        print('\nNo eligió una variable válida')

    return


# Cálculo de a y b
def ab(TC, PC):
    R = 8.314472  # J /molˆ™ K

    a = round(0.45724 * ((R * TC) ** 2 / PC), 5)
    b = 0.07780 * ((R * TC) / PC)

    return round(a, 5), b


# Cálculo de k
def kappa(W):
    return round(0.37464 + (1.54226 * W) - (0.26992 * W ** 2), 5)


# Cálculo de alfa
def alfa(k, TC, T):
    Tr = T/TC
    return round((1 + (k * (1 - Tr ** 0.5))) ** 2, 5)


# Reglas de Mezclado
def abmix(TC, PC, X):
    A = []
    B = []

    for i in range(len(TC)):
        a, b = ab(TC[i], PC[i])
        A.append(a), B.append(b)

    amix = 0
    bmix = 0

    for i in range(len(X)):
        for j in range(len(X)):
            amix = amix + X[i] * X[j] * (A[i] * A[j]) ** 0.5
        bmix = bmix + (X[i] * B[i])

    return round(amix, 5), bmix, A, B


def alfamix(TC, W, X, T):
    ALFA = []
    mix = 0

    for i in range(len(TC)):
        r = alfa(kappa(W[i]), TC[i], T)
        ALFA.append(r)

    for i in range(len(X)):
        for j in range(len(X)):
            mix = mix + X[i] * X[j] * (ALFA[i] * ALFA[j]) ** 0.5

    return round(mix, 5), ALFA


# Presión
def press(tc, pc, w, x, esp, T, v):
    from tabulate import tabulate

    R = 8.314472  # J /molˆ™ K

    a, b, A, B = abmix(tc, pc, x)
    ALFA, M_ALFA = alfamix(tc, w, x, T)

    res = [['Especie', 'Fórmula', 'X', 'TC [K]', 'PC [bar]', 'a', 'b', 'alfa']]

    for i in range(len(A)):
        sub = []
        sub.append(i + 1), sub.append(esp[i]), sub.append(x[i]), sub.append(tc[i]), sub.append(round(pc[i]/1e5,2)), \
            sub.append(A[i]), sub.append(format(B[i], '0.5e')), sub.append(M_ALFA[i])
        res.append(sub)

    print(tabulate(res))

    impresion(a, b, ALFA, 1)

    P = round(((R * T) / (v - b)) - (a * ALFA / ((v * (v + b)) + (b * (v - b)))), 2)
    print(f'Presión: {P} Pa   =  {round(P / 1e5, 3)} bar ')

    M_A_ALFA = [A[i]*M_ALFA[i] for i in range(len(A))]

    cf, f = fugacidad(P, T, v, x, a*ALFA, b, M_A_ALFA, B)

    res = [['Especie', 'Fórmula', 'Coeficiente de Fugacidad', 'Fugacidad [bar]']]

    for i in range(len(A)):
        sub = []
        sub.append(i + 1), sub.append(esp[i]), sub.append(cf[i]), sub.append(f[i])
        res.append(sub)

    print(tabulate(res))

    return


# Temperatura
def temp(tc, pc, w, P, v):
    from scipy.optimize import newton

    R = 8.314472  # J /molˆ™ K

    a, b = ab(tc[0], pc[0])
    k = kappa(w[0])
    impresion(a, b, k, 2)

    g1 = R / (v - b)
    g2 = a / ((v * (v + b)) + (b * (v - b)))

    def f(t):
        return P - (g1 * t) + (g2 * (1 + (k * (1 - (t / tc[0]) ** 0.5))) ** 2)

    T = tc[0]
    T0 = newton(f, T, fprime=None, args=(), tol=1.5e-08, maxiter=50, fprime2=None)
    T0 = round(T0.real, 3)
    print(f'Temperatura (Kelvin): {T0}')

    return


# Volumen
def vol(tc, pc, w, x, esp, P, T):
    import numpy as np
    from tabulate import tabulate

    R = 8.314472  # J /molˆ™ K

    a, b, A, B = abmix(tc, pc, x)
    ALFA, M_ALFA = alfamix(tc, w, x, T)

    res = [['Especie', 'Fórmula', 'X', 'TC [K]', 'PC [bar]', 'a', 'b', 'alfa']]

    for i in range(len(A)):
        sub = []
        sub.append(i + 1), sub.append(esp[i]), sub.append(x[i]), sub.append(tc[i]), sub.append(round(pc[i]/1e5,2)), \
            sub.append(A[i]), sub.append(format(B[i], '0.5e')), sub.append(M_ALFA[i])
        res.append(sub)

    print(tabulate(res))

    impresion(a, b, ALFA, 1)

    n1 = R * T
    n2 = a * ALFA

    c1 = P
    c2 = (P * b) - n1
    c3 = n2 - (3 * P * b ** 2) - (2 * n1 * b)
    c4 = (P * b ** 3) + (n1 * b ** 2) - (n2 * b)
    coeff = [c1, c2, c3, c4]

    v = np.roots(coeff)
    v = v_real(v)
    print(f'Volumen (m^3/mol): {round(v,8)}')

    M_A_ALFA = [A[i] * M_ALFA[i] for i in range(len(A))]

    cf, f = fugacidad(P, T, v, x, a * ALFA, b, M_A_ALFA, B)

    res = [['Especie', 'Fórmula', 'Coeficiente de Fugacidad', 'Fugacidad [bar]']]

    for i in range(len(A)):
        sub = []
        sub.append(i + 1), sub.append(esp[i]), sub.append(cf[i]), sub.append(f[i])
        res.append(sub)

    print(tabulate(res))

    return


# Impresiones
def impresion(a, b, extra, op):
    print('\nParámetros a & b de la ecuación de Peng-Robinson')
    print(' a = ', a, ' [(J m^3) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] ')
    if op == 1:
        print(' alfa =', extra)
    else:
        print(' k =', extra)
    return


# Volumen real
def v_real(vs):
    import numpy as np

    v = []

    for n in vs:
        if np.isreal(n) == True and n > 0:
            v.append(np.real(n))

    return max(v)


# Fugacidad
def fugacidad(p, t, v, x, a_alfa, b, A_ALFA, B):
    import numpy as np

    R = 8.314472  # J /molˆ™ K
    z = p*v/(R*t)

    CF = []
    F = []
    A_ALFA1 = A_ALFA.copy()
    x1 = x.copy()

    for j in range(len(B)):
        mix = 0

        for i in range(len(B)):
            mix += x[i] * (A_ALFA1[j] * A_ALFA[i]) ** 0.5

        bs = B[j]/b

        c1 = bs*(z-1)
        c2 = - np.log((v - b) * p / (R * t))

        c3_1 = a_alfa / (2 ** 1.5 * b * R * t)
        c3_2 = bs - 2 * mix / (a_alfa)
        c3_3 = np.log((v + b*(1 + 2**.5)) / (v + b*(1 - 2**.5)))

        c3 = c3_1 * c3_2 * c3_3

        cf = np.exp(c1 + c2 + c3)
        f = cf * p * x1[j] / 1e5

        A_ALFA.append(A_ALFA[0]), x.append(x[0])
        del A_ALFA[0], x[0]

        CF.append(round(cf, 5)), F.append(round(f, 2))

    return CF, F


# Gráficas: Fugacidad Mezcla Binaria
def fug_bin(TC, PC, W, ESP, T):
    import numpy as np
    from matplotlib import pyplot as plt

    R = 8.314472  # J /molˆ™ K
    P = np.linspace(1, 800, 100) * 1e5
    X1 = [0.001, 0.2, 0.4, 0.6, 0.8, 1]

    CF = []
    CF1 = []
    CF2 = []

    for i in range(len(X1)):
        X = [X1[i], 1 - X1[i]]
        CF_sub1 = []
        CF_sub2 = []

        for j in range(len(P)):
            a, b, A, B = abmix(TC, PC, X)
            ALFA, M_ALFA = alfamix(TC, W, X, T)

            n1 = R * T
            n2 = a * ALFA

            c1 = P[j]
            c2 = (P[j] * b) - n1
            c3 = n2 - (3 * P[j] * b ** 2) - (2 * n1 * b)
            c4 = (P[j] * b ** 3) + (n1 * b ** 2) - (n2 * b)
            coeff = [c1, c2, c3, c4]

            v = np.roots(coeff)
            v = v_real(v)

            M_A_ALFA = [A[r] * M_ALFA[r] for r in range(len(A))]

            cf, f = fugacidad(P[j], T, v, X, a * ALFA, b, M_A_ALFA, B)
            CF_sub1.append(cf[0])
            CF_sub2.append(cf[1])

        CF1.append(CF_sub1)
        CF2.append(CF_sub2)

    CF.append(CF1), CF.append(CF2)

    for j in range(2):
        plt.figure(j + 1)
        leyenda = []

        for i in range(len(X1)):
            plt.plot(P / 1e5, CF[j][i])
            leyenda.append(f'$y_1 = {X1[i]}$')

        if j == 1:
            leyenda.reverse()

        plt.xlabel('Presión [bar]')
        plt.ylabel(f'$φ_{ESP[0]}$')
        plt.suptitle(f'Mezcla {ESP[0]} - {ESP[1]}')
        plt.title(f'T = {T}K')
        plt.legend(leyenda)
        plt.xlim([0, P[-1] / 1e5])
        plt.grid()
        # plt.show()

        ESP.append(ESP[0])
        del ESP[0]

    return
