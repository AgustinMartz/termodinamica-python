# Termodinámica Aplicada con Python

Este repositorio contiene una colección de programas en Python desarrollados como parte de un proyecto educativo en el marco del **Programa Delfín**, dentro del **Verano de la Investigación Científica y Tecnológica del Pacífico**. Los programas permiten resolver problemas clásicos de termodinámica mediante el uso de ecuaciones de estado como:

- Van der Waals
- Redlich-Kwong
- Peng-Robinson

Además, se incluyen versiones para sustancias puras y mezclas, así como integración con bibliotecas externas como **CoolProp** para la obtención automática de propiedades termodinámicas.

---

## 🎯 Objetivo del proyecto

Fomentar el aprendizaje activo de la termodinámica en estudiantes de ingeniería, a través de la programación y resolución numérica de problemas reales. Los programas aquí incluidos fueron desarrollados por estudiantes de licenciatura como parte de sus estancias de investigación, guiados por el Dr. Agustín Martínez Ruvalcaba (Universidad de Guadalajara).

---

## 👥 Créditos

Este proyecto fue desarrollado en el marco del Programa Delfín (2023–2024), durante estancias de investigación en el Departamento de Ingeniería Química de la Universidad de Guadalajara.

### Estudiantes participantes:

- **Nombre Apellido 1** – Universidad Tecnológica de X
- **Nombre Apellido 2** – Universidad Y
- ...
  
### Asesor responsable:
- **Dr. Agustín Martínez Ruvalcaba**  
  Universidad de Guadalajara – CUCEI  
  [agustin.martinez@academicos.udg.mx]  

---

## 📁 Estructura del repositorio

termodinamica-python/
│
├── README.md                 ← Descripción general del proyecto
├── LICENSE                   ← Licencia (ej. MIT)
├── requirements.txt          ← Librerías necesarias
│
├── van_der_waals/            ← Carpeta 1: Ejemplos con Van der Waals
│   ├── vdw_pura.py
│   ├── vdw_mezcla.py (si lo tienes)
│   └── ejemplo_vdw.ipynb (si usas notebooks)
│
├── redlich_kwong/           ← Carpeta 2: Redlich-Kwong (después)
│   ├── rk_pura.py
│   └── rk_mezcla.py
│
├── peng_robinson/           ← Carpeta 3: Peng-Robinson (después)
│   └── ...
│
├── coolprop_integration/    ← Carpeta 4: Integración con CoolProp
│   └── ...
│
└── CREDITOS.md              ← Archivo opcional para reconocer a los autores
