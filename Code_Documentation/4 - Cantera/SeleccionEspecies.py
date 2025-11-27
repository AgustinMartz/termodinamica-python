def especies_manual():
    # Librerías
    import sys

    # Cantidad de especies en la mezcla
    num = int(input('¿De cuántas especies consta la mezcla?'))
    if num <= 0:
        num = 1

    TC = []
    PC = []
    W = []
    X = []
    ESP = []

    for i in range(num):
        esp = input(f'Fórmula o nombre de la especie {i}: ')
        tc = float(input('Temperatura crítica (Kelvin): '))
        pc = float(input('Presión crítica (bar): ')) * 1e5
        w = float(input('Factor acéntrico: '))

        if num != 1:
            X.append(float(input('Escriba la fracción molar de dicho componente: ')))
        else:
            X.append(1)

        TC.append(tc), PC.append(pc), W.append(w), ESP.append(esp.upper())

        error = 0.005
        if abs(sum(X) - 1) > error:
            print('\nLas fracciones molares no suman 1, por favor reinicie')
            sys.exit()

    return TC, PC, W, X, ESP


def especies():
    # Librerías
    import sys

    # Cantidad de especies en la mezcla
    num = int(input('¿De cuántas especies consta la mezcla?'))
    if num <= 0:
        num = 1

    TC = []
    PC = []
    W = []
    X = []
    ESP = []

    for i in range(num):
        # Selección de especies
        print(f'\n------SELECCIÓN DE SUSTANCIA {i + 1}------')
        esp = int(input('1) Agua \n2) Nitrógeno \n3) Oxígeno \n4) Refrigerante 134a'
                        '\n5) Metano \n6) Hidrógeno \n7) Heptano \n8) Otra \n'))

        # Propiedades Críticas
        if 1 <= esp <= 7:
            tc, pc, w, specie = prop_cantera(esp)
        else:
            sus = input('Escriba la fórmula condensada de la especie deseada: ')
            sus = sus.upper()
            tc, pc, w = prop_chemicals(sus)
            specie = sus

        TC.append(tc), PC.append(pc), W.append(w), ESP.append(specie)

        if num != 1:
            X.append(float(input('Escriba la fracción molar de dicho componente: ')))
        else:
            X.append(1)

    error = 0.005
    if abs(sum(X) - 1) > error:
        print('\nLas fracciones molares no suman 1, por favor reinicie')
        sys.exit()
    return TC, PC, W, X, ESP


# Mezcla Binaria
def especies_binaria():

    TC = []
    PC = []
    W = []
    ESP = []
    num = 2

    for i in range(num):
        # Selección de especies
        print(f'\n------SELECCIÓN DE SUSTANCIA {i + 1}------')
        esp = int(input('1) Agua \n2) Nitrógeno \n3) Oxígeno \n4) Refrigerante 134a'
                        '\n5) Metano \n6) Hidrógeno \n7) Dióxido de Carbono \n8) Heptano \n9) Otra \n'))

        # Propiedades Críticas
        if 1 <= esp <= 8:
            tc, pc, w, specie = prop_cantera(esp)
        else:
            sus = input('Escriba la fórmula condensada de la especie deseada: ')
            tc, pc, w = prop_chemicals(sus)
            specie = sus

        TC.append(tc), PC.append(pc), W.append(w), ESP.append(specie)

    return TC, PC, W, ESP


# Propiedades Cantera
def prop_cantera(esp):
    # Librerías
    import cantera as ct
    import numpy as np

    # Especies disponibles
    especies = {1: [ct.Water(), 'H20'],
                2: [ct.Nitrogen(), 'N2'],
                3: [ct.Oxygen(), 'O2'],
                4: [ct.Hfc134a(), '134a'],
                5: [ct.Methane(), 'CH4'],
                6: [ct.Hydrogen(), 'H2'],
                7: [ct.Heptane(), 'C7H16']
                }

    sus = especies[esp][0]

    # Factor acéntrico
    def factor_acentrico(sustancia):
        t = 0.7 * sustancia.critical_temperature
        sustancia.TQ = t, 1
        p = sustancia.P
        return round(-1 - np.log10(p / sustancia.critical_pressure), 5)

    return sus.critical_temperature, sus.critical_pressure, factor_acentrico(sus), especies[esp][1]


# Propiedades Chemicals
def prop_chemicals(sus):
    # Librerías
    import chemicals as ch

    sus = ch.search_chemical(sus)
    return ch.critical.Tc(sus.CAS), ch.critical.Pc(sus.CAS), ch.acentric.omega(sus.CAS)
