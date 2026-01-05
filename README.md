# Temporal Network, Centrality Measures, and Ordinary Differential Equations

![MATLAB](https://img.shields.io/badge/Language-MATLAB-orange.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

## 📌 Overview

This repository contains the source code for the Master's thesis titled **"Temporal network, centrality measures, and ordinary differential equations"**.

The project implements various algorithms to analyze temporal networks, specifically focusing on receiver centralities using Katz and matrix exponential methods. It explores the intersection of network theory and differential equations to model node importance over time.

## 📂 Project Structure

The core implementation is written in **MATLAB**. All source files are located in the `Thesis_code` directory.

### Prerequisites

For the optimal execution of this code, the following MATLAB library is required:
* **[Chebfun](https://www.chebfun.org/)** – An open-source software system for numerical computing with functions.

## 📜 File Descriptions

The `Thesis_code` directory contains the following scripts and functions:

| File Name | Description |
| :--- | :--- |
| **`temporal_node_evolution.m`** | Calculates receiver centralities (Katz and matrix exponential) and generates summary plots representing the temporal evolution of node importance. |
| **`load_FBsoc_dataset.m`** | Loads the `FBsoc` dataset and constructs adjacency matrices according to the specified time granularity (month, week, or day). |
| **`lagrange_cardinal_basis.m`** | Calculates function handles for the Lagrange Cardinal Basis using Barycentric Interpolation for equidistant points. |
| **`interpolation_solver.m`** | Computes the matrix version of the Time-Ordered Exponential based on Section 3.5, utilizing the theoretical framework established in Chapter 3. |
| **`init.m`** | The main initialization script. It serves as the base for all experiments presented in the thesis, divided into specific sections. |
| **`genCoeffMatrix_SP_interval.m`** | Computes the coefficient matrix of $f(t)\Theta(t-s)$ (where $\Theta$ is the Heaviside theta function) within a basis of orthonormal Legendre polynomials. |
| **`exact_sol.m`** | Computes the exact solution of the Time-Ordered Exponential using the matrix exponential for the piecewise constant case, or `ode45` for the continuous case. |
| **`FBsoc.txt`** | The raw text file containing the temporal edge data for the FBsoc dataset. |

## 🚀 Usage

1.  Ensure you have MATLAB installed.
2.  Install the **Chebfun** library and add it to your MATLAB path.
3.  Clone this repository:
    ```bash
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    ```
4.  Navigate to the `Thesis_code` directory.
5.  Run the `init.m` script to execute the experiments.

## ✍️ Author

* **[Filip Miroslav Kucka]** - *Master's Thesis*
* **[Stefano Pozza, Dr., PH.D.]** - *Supervisor and creator of the genCoeffMatrix_SP_interval.m function*
