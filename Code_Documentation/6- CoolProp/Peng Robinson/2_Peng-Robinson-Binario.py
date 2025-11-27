import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import fsolve
import requests
from PIL import Image
from io import BytesIO

image = str(input("¿Ver la lista de compuestos en la base de datos (s/n)?: "))

def download_and_show_image_from_gdrive(file_id):
    try:
        # Create the download URL
        url = f"https://drive.google.com/uc?export=download&id={file_id}"

        # Download the image
        response = requests.get(url)
        response.raise_for_status()  # Check if the download was successful

        # Convert the response content to an image
        image_data = BytesIO(response.content)
        image = Image.open(image_data)

        # Display the image
        image.show()
    
    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
    except IOError as e:
        print(f"Error opening the image: {e}")

file_id = '1EpP1CEbd3ylr39jsWUGrzeIOlpYBQ7b7'  

if image.lower() == "s":
    download_and_show_image_from_gdrive(file_id)

# Universal gas constant (bar·m^3/mol·K)
R = 8.314462618 * 1e-5

substance_1 = str(input("Ingrese el nombre de la primer sustancia en inglés: \n"))

try: 
    # Critical pressure (bar)
    Pc_1 = CP.PropsSI('PCRIT', substance_1) * 1e-5
    # Critical temperature (K)
    Tc_1 = CP.PropsSI('TCRIT', substance_1)
    # Acentric factor (dimensionless)
    w_1 = CP.PropsSI('acentric', substance_1)
except Exception as e:
    print(f'Error: {e}')
    print('Sustancia no encontrada en la base de datos.')

# bar·m^6/mol^2
a_1 = 0.45724 * (R * Tc_1)**2 / Pc_1
# m^3/mol
b_1 = 0.0778 * R * Tc_1 / Pc_1

print("Presión crítica: ", Pc_1, "bar")
print("Temperatura crítica: ", Tc_1, "K")

#
#
#
#
#

substance_2 = str(input("Ingrese el nombre de la segunda sustancia en inglés: \n"))

try: 
    # Critical pressure (bar)
    Pc_2 = CP.PropsSI('PCRIT', substance_2) * 1e-5
    # Critical temperature (K)
    Tc_2 = CP.PropsSI('TCRIT', substance_2)
    # Acentric factor (dimensionless)
    w_2 = CP.PropsSI('acentric', substance_2)
except Exception as e:
    print(f'Error: {e}')
    print('Sustancia no encontrada en la base de datos.')

# bar·m^6/mol^2
a_2 = 0.45724 * (R * Tc_2)**2 / Pc_2
# m^3/mol
b_2 = 0.0778 * R * Tc_2 / Pc_2

print("Presión crítica: ", Pc_2, "bar")
print("Temperatura crítica: ", Tc_2, "K")

#
#
#
#
#

y1 = float(input(f"Fracción mol de {substance_1}"))
y2 = 1 - y1

amix = y1**2 * a_1 + y2**2 * a_2 + 2 * y1 * y2 * (a_1 * a_2)**0.5*(1-0.092)
bmix = y1 * b_1 + y2 * b_2 
Tpc  = y1 * Tc_1 + y2 * Tc_2
wmix = y1 * w_1 + y2 * w_2
kmix = 0.37464 + 1.54226 * wmix - 0.26992 * wmix**2

print(" amix: ", amix, "bar·m^6/mol^2")
print(" bmix: ", bmix, "m^3/mol")
print(" Temperatura pseudocrítica: ", Tpc, "K")
print(" kmix: ", kmix)


variable = str(input("Variable de estado (P, V o T): \n"))

if variable.lower() == "v":
    
    P = float(input("Presión (bar): "))
    T = float(input("Temperatura (K): "))
    Tr = T / Tpc
    alpha = (1 + kmix*(1 - Tr**0.5))**2
    
    def peng_robinson_v(V):
        term1 = R * T / (V - bmix)
        term2 = (amix * alpha) / (V * (V + bmix) + bmix * (V - bmix))
        return P - (term1 - term2)
    
    V_initial_guess = [bmix + 1e-5, R * T / P, 10 * bmix, 50 * bmix, 100 * bmix]
    
    V_m_roots = []
    for guess in V_initial_guess:
        try:
            root = fsolve(peng_robinson_v, guess)[0]
            if root > 0 and np.isreal(root):
                V_m_roots.append(root)
        except RuntimeError:
            # Handle cases where fsolve does not converge
            pass

    # Filter out duplicates
    V_m_roots = np.unique(np.round(np.array(V_m_roots),6))
    
    print("Volumen(es) en m^3/mol:")
    print(V_m_roots)
    
elif variable.lower() == "p":
    
    V = float(input("Volumen (m^3/mol): "))
    T = float(input("Temperatura (K): "))
    Tr = T / Tpc
    alpha = (1 + kmix*(1 - Tr**0.5))**2
    
    P = R*T/ (V - bmix) + amix*alpha / (V*(V+bmix) + bmix*(V-bmix))
    
    print(P, "bar")
    
elif variable.lower() == "t":
    
    V = float(input("Volumen (m^3/mol): "))
    P = float(input("Presión (bar): "))
    
    
    def peng_robinson_t(T):
        term1 = R * T / (V - bmix)
        term2 = (amix * (1 + kmix*(1 - (T/Tpc)**0.5)**2)) / (V * (V + bmix) + bmix * (V - bmix))
        return P - (term1 - term2)
    
    
    
    T_initial_guess = [1/2 * Tpc, Tpc, 2 * Tpc, 3 * Tpc, 4 * Tpc]
    
    T_m_roots = []
    for guess in T_initial_guess:
        try:
            root = fsolve(peng_robinson_t, guess)[0]
            if root > 0 and np.isreal(root):
                T_m_roots.append(root)
        except RuntimeError:
            # Handle cases where fsolve does not converge
            pass

    # Filter out duplicates
    T_m_roots = np.unique(np.round(np.array(T_m_roots),6))
    
    print("Temperatura(s) en K: ")
    print(T_m_roots)
 
else: 
    print("Variable inválida")



