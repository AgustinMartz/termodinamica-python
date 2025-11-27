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

substance = str(input("Ingrese el nombre de la sustancia en inglés: \n"))

try: 
    # Critical pressure (bar)
    Pc = CP.PropsSI('PCRIT', substance) * 1e-5
    # Critical temperature (K)
    Tc = CP.PropsSI('TCRIT', substance)
    # Acentric factor (dimensionless)
    w = CP.PropsSI('acentric', substance)
    # Universal gas constant (bar·m^3/mol·K)
    R = 8.314462618 * 1e-5
except Exception as e:
    print(f'Error: {e}')
    print('Sustancia no encontrada en la base de datos.')

# bar·m^6/mol^2
a = 0.45724 * (R * Tc)**2 / Pc
# m^3/mol
b = 0.0778 * R * Tc / Pc
# Dimensionless
k = 0.37464 + 1.54226 * w - 0.26992 * w**2

print("Presión crítica: ", Pc, "bar")
print("Temperatura crítica: ", Tc, "K")
print("Factor acentrico: ", w)
print("a: ", a, "bar·m^3/mol*K")
print("b: ", b, "m^3/mol\n")

variable = str(input("Variable de estado (P, V o T): \n"))

if variable.lower() == "v":
    
    P = float(input("Presión (bar): "))
    T = float(input("Temperatura (K): "))
    Tr = T / Tc
    alpha = (1 + k*(1 - Tr**0.5))**2
    
    def peng_robinson_v(V):
        term1 = R * T / (V - b)
        term2 = (a * alpha) / (V * (V + b) + b * (V - b))
        return P - (term1 - term2)
    
    V_initial_guess = [b + 1e-5, R * T / P, 10 * b, 50 * b, 100 * b]
    
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
    Tr = T / Tc
    alpha = (1 + k*(1 - Tr**0.5))**2
    
    P = R*T/ (V - b) + a*alpha / (V*(V+b) + b*(V-b))
    
    print("Presión")
    print(P, "bar")
    
elif variable.lower() == "t":
    
    V = float(input("Volumen (m^3/mol): "))
    P = float(input("Presión (bar): "))
    
    
    def peng_robinson_t(T):
        term1 = R * T / (V - b)
        term2 = (a * (1 + k*(1 - (T/Tc)**0.5)**2)) / (V * (V + b) + b * (V - b))
        return P - (term1 - term2)
    
    
    
    T_initial_guess = [1/2 * Tc, Tc, 2 * Tc, 3 * Tc, 4 * Tc]
    
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
