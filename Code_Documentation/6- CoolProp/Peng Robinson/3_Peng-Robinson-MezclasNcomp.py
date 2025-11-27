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

# Ask the user how many substances they want to input
num_substances = int(input("Número de sustancias: \n"))

# Initialize a list to store information about each substance
substances = []

for i in range(num_substances):
    substance_name = str(input(f"Ingrese el nombre de la sustancia {i+1} en inglés: \n"))
    try: 
        # Get critical properties
        Pc = CP.PropsSI('PCRIT', substance_name) * 1e-5         # Critical pressure in bar
        Tc = CP.PropsSI('TCRIT', substance_name)                # Critical temperature in K
        w = CP.PropsSI('acentric', substance_name)              # Acentric factor (dimensionless)
        
        # Calculate 'a' and 'b' constants
        a = 0.45724 * (R * Tc)**2 / Pc                          # bar·m^6/mol^2
        b = 0.0778 * R * Tc / Pc                                # m^3/mol

        # Store the information in a dictionary
        substance_info = {
            "name": substance_name,
            "Pc": Pc,
            "Tc": Tc,
            "w": w,
            "a": a,
            "b": b
        }
        substances.append(substance_info)

    except Exception as e:
        print(f'Error: {e}')
        print(f'Sustancia {substance_name} no encontrada en la base de datos.')
        
# Ask the user for the mole fractions
mole_fractions = []
for i in range(num_substances):
    x_i = float(input(f"Ingrese la fracción molar de la sustancia {substances[i]['name']}: \n"))
    mole_fractions.append(x_i)

# Calculate mix properties
Tpc  = sum(mole_fractions[i] * substances[i]['Tc'] for i in range(num_substances))
bmix = sum(mole_fractions[i] * substances[i]['b'] for i in range(num_substances))
wmix = sum(mole_fractions[i] * substances[i]['w'] for i in range(num_substances))
kmix = 0.37464 + 1.54226 * wmix - 0.26992 * wmix**2

# amix = sum(mole_fractions[i] * substances[i]['a'] for i in range(num_substances))
amix = 0
for i in range(num_substances):
    for j in range(num_substances):
        aij   = (substances[i]['a'] * substances[j]['a'])**0.5
        amix += mole_fractions[i] * mole_fractions[j] * aij

for substance in substances:
    print("\nInformación de la sustancia:")
    print(f"Nombre: {substance['name']}")
    print(f"Presión crítica: {substance['Pc']} bar")
    print(f"Temperatura crítica: {substance['Tc']} K\n")

print("Propiedades de la mezcla: ")
print(" amix: ", amix, "bar·m^6/mol^2")
print(" bmix: ", bmix, "m^3/mol")
print(" Temperatura pseudocrítica: ", Tpc, "K")
print(" kmix: ", kmix,"\n")
    
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