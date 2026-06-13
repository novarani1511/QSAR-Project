# QSAR and Molecular Modeling Studies on Histamine H1-Receptor Antagonists

This repository contains the complete pipeline for a Quantitative Structure-Activity Relationship (QSAR) and Molecular Docking study on Histamine H1-receptor antagonists. The project was developed as a comprehensive thesis-level study using Python, RDKit, and AutoDock Vina.

## Project Overview

The main objective of this study is to predict the biological activity (pIC50) of Histamine H1 antagonists using Multiple Linear Regression (MLR) and to explore their binding interactions at the atomic level through molecular docking.

### Workflow Pipeline
1. **Data Collection (`01_data_collection.py`)**: Fetches bioactivity data for the human Histamine H1 receptor (CHEMBL231) directly from the ChEMBL API.
2. **Data Curation (`02_data_cleaning.py`)**: Filters invalid data, limits standard units to nanomolar (nM), converts IC50 to pIC50, and handles duplicate SMILES by averaging their activities.
3. **Descriptor Calculation (`03_descriptor_calculation.py`)**: Generates 3D conformers, performs geometry optimization using the MMFF94 force field, and calculates 219 standard RDKit descriptors.
4. **QSAR Modeling (`04_qsar_modeling.py`)**: Removes low-variance and highly correlated descriptors (>0.85). Performs Recursive Feature Elimination (RFE) to select the top 5 descriptors to avoid overfitting in the Multiple Linear Regression (MLR) model.
5. **Applicability Domain (`05_applicability_domain.py`)**: Evaluates the model's reliability using a Williams Plot to detect structural outliers and response outliers.
6. **Molecular Docking (`06_molecular_docking.py`)**: Prepares the most active ligand using `meeko`, downloads the 3RZE crystal structure, and automates docking via AutoDock Vina.

## Installation and Requirements

To run this project locally, ensure you have Python 3.9+ installed.

```bash
pip install pandas numpy scikit-learn matplotlib seaborn rdkit chembl_webresource_client meeko
```

*Note: For the molecular docking step, you must manually download the [AutoDock Vina executable](https://github.com/ccsb-scripps/AutoDock-Vina/releases) and place `vina.exe` in the root folder. Additionally, you will need AutoDockTools (ADT) to manually prepare the receptor (`receptor.pdbqt`) from the raw `3rze.pdb`.*

## Results

**Multiple Linear Regression (MLR) Statistics:**
- **R² Training**: 0.286
- **Q² (5-Fold CV)**: 0.257
- **R² Test (External)**: 0.346

**QSAR Equation:**
```text
pIC50 = 6.575 - 0.496 * [BCUT2D_CHGHI] - 0.332 * [BCUT2D_CHGLO] + 0.313 * [VSA_EState9] + 0.501 * [fr_Ndealkylation2] + 0.495 * [fr_bicyclic]
```

**Applicability Domain:**
The Williams Plot (`williams_plot.png`) confirms the predictive capability boundaries of the model, with the leverage warning limit ($h^*$) calculated at 0.044.

## Purpose
Developed for academic and thesis research purposes to demonstrate a fully reproducible, programmatic workflow in computational chemistry.
