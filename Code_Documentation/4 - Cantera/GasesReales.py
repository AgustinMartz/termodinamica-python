# Librerías
import SeleccionEspecies as SE
import E_VDW as VDW
import E_Redlich as RED
import E_Peng_Robinson as PENG

# Opción única a volúmenes
# Sí - 1 // No - 0
volume = 1

# Entrada de datos manual
# Sí - 1 // No - 0
manual = 0

# Selección de especies
if manual == 0:
    # Auto
    TC, PC, W, X, ESP = SE.especies()
else:
    TC, PC, W, X, ESP = SE.especies_manual()

# Selección de Ecuación
print('\n------SELECCIÓN DE ECUACIÓN------')
ec = int(input('1) Van der Waals \n2) Redlich-Kwong \n3) Peng-Robinson \n'))
print()

if ec == 1:
    VDW.vdw(PC, TC, X, volume, ESP)
elif ec == 2:
    RED.redlich(PC, TC, X, volume, ESP)
elif ec == 3:
    PENG.peng_rob(PC, TC, W, X, volume, ESP)
else:
    print('Fuera de rango, por defecto se utilizará Peng-Robinson\n')
    PENG.peng_rob(PC, TC, W, X, volume, ESP)
