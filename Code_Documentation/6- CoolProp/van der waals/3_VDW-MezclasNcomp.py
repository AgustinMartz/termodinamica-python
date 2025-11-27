import CoolProp.CoolProp as CP
import numpy as np

# Constante universal de los gases en J/(mol*K)
R = 8.3144621

# Definimos las constantes
def van_der_waals_constants(Tc, Pc, x, k):
    n = len(Tc)
    a = np.zeros((n, n))
    b = np.zeros(n)

# Calculamos las constantes y se repite para cada uno
    for i in range(n):
        b[i] = (R * Tc[i]) / (8 * Pc[i])
        for j in range(n):
            a[i, j] = np.sqrt((27 / 64) * (R ** 2) * (Tc[i] * Tc[j]) / (Pc[i] * Pc[j])) * (1 - k[i, j])

    a_m = np.sum([x[i] * x[j] * a[i, j] for i in range(n) for j in range(n)])
    b_m = np.sum([x[i] * b[i] for i in range(n)])

    return a_m, b_m


print('ECUACIÓN DE ESTADO DE VAN DER WAALS PARA MEZCLAS MULTICOMPONENTE')

# Obtener el número de componentes
n = int(input('Ingrese el número de componentes en la mezcla: '))

# Seleccionar las sustancias y fracciones molares
sustancias = []
x = []
for i in range(n):
    sustancia = input(f'Escribe el nombre de la sustancia {i + 1} en inglés: ')
    sustancias.append(sustancia)
    fraccion_molar = float(input(f'Ingrese la fracción molar de la sustancia {i + 1} (entre 0 y 1): '))
    x.append(fraccion_molar)

# Verificar que la suma de las fracciones molares sea 1
if not np.isclose(sum(x), 1.0):
    raise ValueError('La suma de las fracciones molares debe ser igual a 1.')

# Obtener las propiedades críticas
Tc = []
Pc = []
try:
    for sustancia in sustancias:
        Tc.append(CP.PropsSI('Tcrit', sustancia))
        Pc.append(CP.PropsSI('Pcrit', sustancia))

    # Coeficientes de interacción binaria (suponiendo todos ceros para simplificación)
    k = np.zeros((n, n))

    # Calcular las constantes de Van der Waals para la mezcla
    a_m, b_m = van_der_waals_constants(Tc, Pc, x, k)

    print('a_m = ', a_m, ' [(J m^3) / mol^2] ')
    print('b_m = ', format(b_m, '0.5e'), ' [m^3 / mol] \n')

    # Seleccionar la variable a calcular
    tipo = input('''Elija la variable que desea calcular :
    Teclear P si va a calcular la presión
    Teclear T si va a calcular la temperatura
    Teclear V si va a calcular el volumen específico \n''')

    if tipo.lower() == 'p':
        T = float(input('Temperatura (en Kelvin): '))
        V = float(input('Volumen específico (en m^3/mol): '))
        P = ((R * T) / (V - b_m)) - (a_m / V ** 2)
        P = round(P, 2)
        print('La presión es igual a ', P, ' Pa   =  ', round(P / 100000, 3), ' bar ')

    elif tipo.lower() == 't':
        P = float(input('Presión (en bar): '))
        P = P * 100000  # Convertir bar a Pa
        V = float(input('Volumen específico (en m^3/mol): '))
        T = (P + (a_m / V ** 2)) * ((V - b_m) / R)
        T = round(T, 2)
        print('La temperatura es igual a ', T, ' K ')

    elif tipo.lower() == 'v':
        P = float(input('Presión (en bar): '))
        P = P * 100000  # Convertir bar a Pa
        T = float(input('Temperatura (en Kelvin): '))

        # Resolver la ecuación cúbica para V usando numpy.roots
        coeff = [P, -(R * T + P * b_m), a_m, -a_m * b_m]
        V_roots = np.roots(coeff)
        V_roots = np.real(V_roots[np.isreal(V_roots)])  # Filtrar raíces reales

        print('Volúmenes específicos posibles:')
        for i, v in enumerate(V_roots, start=1):
            print(f'v{i} = {np.round(v, 6)} [m^3/mol]')

    else:
        print('No eligió una variable válida')
except Exception as e:
    print(f'Error al obtener las propiedades de la sustancia: {e}')
