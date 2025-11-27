# Principal
def redlich(pc, tc, x, volumen, esp):
    from tabulate import tabulate

    print('------ECUACIÓN DE ESTADO DE  REDLICH-KWONG------')

    # Parámetros de mezclado amix y bmix
    a, b, A, B = abmix(tc, pc, x)

    res = [['Especie', 'Fórmula', 'X', 'TC [K]', 'PC [bar]', 'a', 'b']]

    for i in range(len(A)):
        sub = []
        sub.append(i + 1), sub.append(esp[i]), sub.append(x[i]), sub.append(tc[i]), sub.append(round(pc[i]/1e5,2)), \
            sub.append(A[i]), sub.append(format(B[i], '0.5e'))
        res.append(sub)

    print(tabulate(res))

    print('\nParámetros generales a & b de la ecuación de Redlich-Kwong')
    print(' a = ', a, ' [(J m^3 K^0.5) / mol^2] ')
    print(' b = ', format(b, '0.5e'), ' [m^3 / mol] \n ')

    if volumen == 1:
        tipo = 'v'
    else:
        tipo = input('''\nElija la varible que desea calcular :
            Teclear P si va a calcular la presión
            Teclear T si va a calcular la temperatura
            Teclear v si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        v = float(input('Volumen específico (en m^3/mol): '))
        P = press(a, b, T, v)

        print(f'Presión: {P} Pa   =  {round(P / 1e5, 3)} bar ')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): ')) * 1e5
        v = float(input('Volumen específico (en m^3/mol): '))
        T = temp(a, b, tc[0], P, v)

        print(f'Temperatura (Kelvin): {T}')

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): ')) * 1e5
        T = float(input('Temperatura (en Kelvin): '))
        v = vol(a, b, P, T)

        print(f'Volumen (m^3/mol): {round(v, 8)}')

    else:
        print('\nNo eligió una variable válida')

    cf, f = fugacidad(P, T, v, x, a, b, A, B)

    if len(x) == 1:
        fug_sp(tc, pc, esp, x)
    elif len(x) == 2:
        fug_bin(tc, pc, esp, T)

    res = [['Especie', 'Fórmula', 'Coeficiente de Fugacidad', 'Fugacidad [bar]']]

    for i in range(len(A)):
        sub = []
        sub.append(i + 1), sub.append(esp[i]), sub.append(cf[i]), sub.append(f[i])
        res.append(sub)

    print(tabulate(res))

    return

# Cálculo de a y b
def ab(TC, PC):
    R = 8.314472  # J /molˆ™ K

    a = (9 * ((2 ** (1 / 3)) - 1)) ** (-1) * ((R ** 2 * TC ** 2.5) / PC)
    b = (((2 ** (1 / 3)) - 1) / 3) * ((R * TC) / PC)

    return round(a, 5), b


# Regla de Mezclado
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


# Presión
def press(a, b, T, v):
    R = 8.314472  # J /molˆ™ K
    P = round(((R * T) / (v - b)) - (a / (T ** 0.5 * v * (v + b))), 2)

    return P


# Temperatura
def temp(a, b, tc, P, v):
    from scipy.optimize import newton

    R = 8.314472  # J /molˆ™ K

    m1 = -(P * (v - b)) / R
    m2 = -a / (R * v)

    def f(t):
        return t ** (3 / 2) + (m1 * t ** 0.5) + m2

    T = tc
    T0 = newton(f, T, fprime=None, args=(), tol=1.5e-08, maxiter=80, fprime2=None)
    T0 = round(T0.real, 2)

    return T0


# Volumen
def vol(a, b, P, T):
    import numpy as np

    R = 8.314472  # J /molˆ™ K

    c1 = 1
    c2 = -(R * T) / P
    c3 = (a / (P * T ** 0.5)) - (b ** 2) - ((b * R * T) / P)
    c4 = -(a * b) / (P * T ** 0.5)
    coeff = [c1, c2, c3, c4]

    v = np.roots(coeff)
    v = v_real(v)

    return v


# Volumen real
def v_real(vs):
    import numpy as np

    v = []

    for n in vs:
        if np.isreal(n) == True and n > 0:
            v.append(np.real(n))

    return max(v)


# Fugacidad
def fugacidad(p, t, v, x, a, b, A, B):
    import numpy as np

    R = 8.314472  # J /molˆ™ K
    z = p*v/(R*t)

    CF = []
    F = []
    A1 = A.copy()
    x1 = x.copy()

    for j in range(len(A)):
        mix = 0

        for i in range(len(A)):
            mix += x[i] * (A1[j] * A[i]) ** 0.5

        bs = B[j]/b

        c1 = bs*(z-1)
        c2 = - np.log((v - b) * p / (R * t))

        c3_1 = np.log(1 + b/v) / (b * R * t ** 1.5)
        c3_2 = bs*a - 2*mix

        c3 = c3_1 * c3_2

        cf = np.exp(c1 + c2 + c3)
        f = cf * p * x1[j] / 1e5

        A.append(A[0]), x.append(x[0])
        del A[0], x[0]

        CF.append(round(cf, 5)), F.append(round(f, 2))

    return CF, F


# Gráficas: Fugacidad Mezcla Binaria
def fug_bin(TC, PC, ESP, T):
    import numpy as np
    from matplotlib import pyplot as plt

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

            v = vol(a, b, P[j], T)

            cf, f = fugacidad(P[j], T, v, X, a, b, A, B)
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


# Gráficos: Fugacidad especie pura
def fug_sp(TC, PC, ESP, X):
    import numpy as np
    from matplotlib import pyplot as plt

    P = np.linspace(1, 800, 100) * 1e5
    T = np.linspace(300, 700, 5)

    CF = []

    for i in T:
        CF1 = []

        for j in P:
            a, b, A, B = abmix(TC, PC, X)

            v = vol(a, b, j, i)

            cf, f = fugacidad(j, i, v, X, a, b, A, B)
            CF1.append(cf[0])

        CF.append(CF1)

    leyenda = []

    for i in range(len(T)):
        plt.plot(P / 1e5, CF[i])
        leyenda.append(f'$T = {T[i]}K$')

    plt.xlabel('Presión [bar]')
    plt.ylabel(f'$φ_{ESP[0]}$')
    plt.title(f'Fugacidad de {ESP[0]}')
    plt.legend(leyenda)
    plt.xlim([0, P[-1] / 1e5])
    plt.grid()
    # plt.show()

    return
