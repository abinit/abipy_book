# Lumabi: Luminescence Lineshape of Point Defects

This collection of chapters introduces **Lumabi**, a collection of Python modules integrated in the [AbiPy](https://github.com/abinit/abipy) framework designed to compute the luminescence lineshape of point defects in solids using [Abinit](https://www.abinit.org/).

## Overview

Lumabi is structured around four main Python modules, each automating a key stage of the computational workflow for phonon-resolved luminescence spectra of defects. These modules are designed to be used together, the output of one module serving as input for another, or separately.

### 1. LumiWork Module

![The LumiWork module, an AbiPy Workflow that automates ABINIT DFT tasks with $\Delta$SCF constrained occupations.](../paper/LumiWork.pdf)

The LumiWork module automates ABINIT DFT tasks with $\Delta$SCF constrained occupations. It manages two structural relaxations (ground and excited states), and optional static SCF and non-SCF band structure calculations. The main outputs are collected in netcdf format, ready for post-processing.

### 2. $\Delta$SCF Post-Processing Module

![The $\Delta$SCF module, designed to post-process $\Delta$SCF constrained-occupation calculations using a one-dimensional configuration-coordinate model.](../paper/dSCF_post_process.pdf)

This module processes the netcdf output files from LumiWork, analyzing them with a one-dimensional configuration coordinate model (1D-CCM). It computes transition energies, Huang-Rhys factors, effective phonon frequencies, lineshapes using the 1D-CCM, and helps analyze atomic relaxations.

### 3. IFCs Embedding Module

![The IFCs embedding module, allowing to calculate defect phonons in large supercells.](../paper/IFCs_embedding.pdf)

The IFCs Embedding module enables the calculation of defect phonons in large supercells by combining interatomic force constants (IFCs) from both pristine and defect systems. This approach allows one to captures both the coupling with long-wavelength and localized phonon modes.

### 4. Lineshape Calculation Module

![The lineshape module, allowing to compute the temperature-dependent spectra.](../paper/lineshape.pdf)

The Lineshape module computes the Huang-Rhys spectral function and generates (temperature-dependent) photoluminescence spectra using the generating function approach. It takes as input the zero-phonon line energy, atomic displacements or forces, and phonon modes (potentially from the IFCs embedding module), and produces the final luminescence spectrum.


## Tutorial Structure
This tutorial is organized as follows:

1. **Theory**: Presents the formalism for luminescence lineshape calculations, from Fermi's golden rule to practical models for solids, and introduces the Huang-Rhys theory and generating function approach.
2. **LumiWork Workflow**: Guides through setting up and running the automated workflow for $\Delta$SCF calculations.
3. **1D Post-Processing**: Shows how to analyze results using the single effective phonon mode, also called one-dimensional configuration coordinate model (1D-CCM).
4. **Multi-Phonon Lineshape**: Explains how to compute the full phonon-resolved luminescence spectrum using the generating function approach, and how to analyze the localization of phonon modes.
5. **IFC Embedding**: Details the embedding approach for obtaining phonon modes in large supercells with defects.


