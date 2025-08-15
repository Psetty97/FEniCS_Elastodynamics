# FEniCS_Elastodynamics

Implementation of **explicit** and **implicit** *Asynchronous Variational Integrators (AVI)* for solving **2D elastodynamics** problems using **Finite Element Method (FEM)** in [FEniCS](https://fenicsproject.org/).

---

## 📌 Overview
This repository contains Python scripts and mesh files for simulating 2D elastodynamics using both **explicit** and **implicit AVI schemes**.  
The methods are formulated based on the **variational principles of mechanics** and are designed to handle problems with **multiple time scales** efficiently, utilising a priority queue implementation in Python.

---

## 📂 Repository Structure
FEniCS_Elastodynamics/
- │
- ├── AVI_Implicit/                      # Scripts for implicit AVI method
- │   └── (files implementing implicit integration)
- │
- ├── AVI_Elastodynamics_explicit.py     # Explicit AVI implementation
- ├── square_lamina_0.1.msh               # Gmsh mesh file
- ├── square_lamina_0.1.xml               # FEniCS XML mesh file
- └── README.md                           # Project documentation

---

## ⚙️ Features
- **Explicit AVI implementation** for 2D elastodynamics
- **Implicit AVI implementation** for increased stability
- Mesh generation using Gmsh (`.msh` and `.xml` formats)
- Energy-preserving integration for long-term stability

---

## 📊 Problem Setup
- **Domain:** Square lamina mesh (`square_lamina_0.1.msh`)
- **Material model:** Hyper elasticity Neo-Hookean model
- **Boundary conditions:** Dirichlet
- **Time integration:** Explicit & implicit AVI

---
### 1️⃣ Prerequisites
Install **FEniCS** (tested on version ≥ 2019.1.0):

### 📈 References
- Lew, A., Marsden, J. E., Ortiz, M., & West, M. (2003). Asynchronous variational integrators. DOI 10.1007/s00205-002-0212-y.
