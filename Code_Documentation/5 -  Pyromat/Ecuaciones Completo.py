from Ecuaciones_Función import *

modo = input('Bienvenido, deseas trabajar en modo: \n'
             '1)Automático \n'
             '2)Manual \n'
             '')
sistema = input('Deseas tabajar con: \n'
                '1) Sustancias puras \n'
                '2) Mezclas \n'
                ':')
if sistema == "2":
    opcion = input('Ingresa la ecuación que deseas utilizar: \n'
                   '1) Van der Waals \n'
                   '2) Redlich-Kwong \n'
                   '3) PengRobinson \n'
                   ':')
elif sistema=="1":
    ecuacion = input("Bienvenido, por favor introduce la ecuación de estado que deseas utilizar: \n"
                     "1) Van der Waals \n"
                     "2) Redlich-Kwong \n"
                     "3) Peng-Robinson \n"
                     "4) Dieterichi \n"
                     "5) Lee_Kesler \n"
                     ":")
if modo == "1":
    if sistema == "1":
        if ecuacion == "1":
            vanwaals()

        elif ecuacion == "2":
            redlich()

        elif ecuacion == "3":
            pengro()

        elif ecuacion == "4":
            dieteri()

        elif ecuacion == "5":
            lekes()

    elif sistema == "2":
        if opcion == "1":
            vandermezcla()
        elif opcion == "2":
            redlich_mezcla()
        elif opcion == "3":
            pengMezcla()

elif modo == "2":
    if sistema == "1":
        if ecuacion == "1":
            vanwaalsmanual()
        elif ecuacion == "2":
            redlichmanual()
        elif ecuacion == "3":
            pengromanual()
        elif ecuacion == "4":
            lekesmanual()
        elif ecuacion == "5":
            dieterimanual()

    elif sistema == "2":
        if opcion == "1":
            vanwaalsmanual()
        elif opcion == "2":
            redlich_mezclamanual()
        elif opcion == "3":
            pengMezclaManual()












