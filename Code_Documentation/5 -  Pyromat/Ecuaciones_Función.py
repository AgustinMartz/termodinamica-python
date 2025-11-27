import chemicals as cm
import pyromat as pm
import math as math
from VanderWaals import vanwaals
from VanderManual import vanwaalsmanual
from RedlichManual import redlichmanual
from Pengmanual import pengromanual
from LeeKManual import lekesmanual
from DieterichiManual import dieterimanual
from RedlichKwong import redlich
from PengRobinson import pengro
from Dieterichi import dieteri
from Lee_Kesler import lekes
from VanderWaals_mezclas import vandermezcla
from RedlichMezcla import redlich_mezcla
from PengRobinsonMezcla import pengMezcla
from VanderManual import vanwaalsmanual
from RedlichMezclaManual import redlich_mezclamanual
from PengroMezclaM import  pengMezclaManual

def ecuacionesfun():
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
    elif sistema == "1":
        ecuacion = input("Bienvenido, por favor introduce la ecuación de estado que deseas utilizar: \n"
                         "1) Van der Waals \n"
                         "2) Redlich-Kwong \n"
                         "3) Peng-Robinson \n"
                         "4) Dieterichi \n"
                         "5) Lee_Kesler \n"
                         ":")

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






