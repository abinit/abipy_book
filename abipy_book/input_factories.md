---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Factory functions

Abipy provides factory functions to build input files for typical calculations.
These functions return `AbinitInput` or `MultiDataset` objects, depending
on the number of steps required by the calculation.

One can use the factories to generate automatically input files or
call these functions inside python code to build workflows for high-throughput applications.
Note that the default values do not always correspond to the default behaviour of Abinit.
In particular, the majority of the factory functions construct input files 
for **spin-polarized calculations** (`nsppol=2`) with a **Fermi-Dirac** occupation scheme and 
a physical temperature of **0.1 eV**. 
It is always possible to change the default behaviour either
by passing these options to the factory function or by changing the objects returned by the factory.

Also note that the factory functions do not use `get*` or `ird*` variables to connect the different 
steps. Client code is in charge of connecting the different parts.
For a command line interface, use the `abinp.py` script.

```{code-cell} ipython3
import os
import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

import abipy.data as abidata
import abipy.abilab as abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook
from abipy.abilab import AbinitInput
```

### Ground-state calculation

Let us generate an input file for a standard GS calculation for silicon in which 
the structure is read from an external CIF file:

```{code-cell} ipython3
si_cif = abidata.cif_file("si.cif")
pseudos = os.path.join(abidata.pseudo_dir, "14si.pspnc")

# Build input for GS calculation (unpolarized, no smearing, 1000 k-points per reciprocal atom) 
# ecut must be specified because this pseudopotential does not provide hints for ecut.
# kppa stands for k-point per reciprocal atom.
gs_inp = abilab.gs_input(
    si_cif, pseudos,
    kppa=1000, ecut=8, spin_mode="unpolarized", smearing=None) # change default

gs_inp.set_mnemonics(True)
gs_inp
```

### Input variables for band structure calculation + DOS

A slightly more complicated example:

```{code-cell} ipython3
---
code_folding: []
run_control:
  marked: true
---
# GS run + NSCF on a path + NSCF run on a k-mesh to compute the DOS
multi = abilab.ebands_input(si_cif, pseudos,
                            ecut=8, spin_mode="unpolarized", smearing=None, dos_kppa=5000)

multi
```

### Factories for GW calculations

```{code-cell} ipython3
# Generate an input file for GW calculations with the plasmon-pole model.
# The calculations consists of a GS run to get the density followed by a 
# nscf-run to compute the WFK file with `nscf_nband` states.
# The cutoff for the screening is given by `ecuteps` while the cutoff for
# the exchange part of the self-energy is equal to ecut.
# kppa defines the k-point sampling.
kppa = 1000
ecut = ecutsigx = 8
ecuteps = 2
nscf_nband = 50

multi = abilab.g0w0_with_ppmodel_inputs(
    si_cif, pseudos, kppa, nscf_nband, ecuteps, ecutsigx,
    ecut=ecut, smearing=None, spin_mode="unpolarized")

multi
```
