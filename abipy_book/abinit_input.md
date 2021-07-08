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

# The AbinitInput object

The creation of the Abinit input file is one of the most repetive and error-prone operations
we have to perform before running our calculations.
To facilitate the creation of the input files, AbiPy provides the `AbinitInput` object,
a dict-like object storing the Abinit variables and providing methods to automate
the specification of multiple parameters.

This notebook discusses how to create an `AbinitInput` and how to define the parameters of the calculation.
In the last part, we introduce the `MultiDataset` object that is mainly designed for the generation
of multiple inputs sharing the same structure and the same list of pseudopotentials.
In another [notebook](./input_factories.ipynb), we briefly discuss how to use factory functions
to generate automatically input objects for typical calculations.

## Creating an AbinitInput object

```{code-cell} ipython3

import os
import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

import abipy.data as abidata
import abipy.abilab as abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook
from abipy.abilab import AbinitInput

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

To create an Abinit input, we must specify the paths of the pseudopotential files.
In this case, we have a pseudo named `14si.pspnc` located
in the `abidata.pseudo_dir` directory.

```{code-cell} ipython3
inp = AbinitInput(structure=abidata.cif_file("si.cif"),
                  pseudos="14si.pspnc", pseudo_dir=abidata.pseudo_dir)
```

`print(inp)` returns a string with our input.
In this case, the input is almost empty since only the structure and the pseudos have been specified.

```{code-cell} ipython3
print(inp)
```

Inside the jupyter notebook, it is possible to visualize the input in HTML including the links to the official ABINIT documentation:

```{code-cell} ipython3
inp
```

The input *has* a structure:

```{code-cell} ipython3
print(inp.structure)
```

and a list of `Pseudo` objects:

```{code-cell} ipython3
for pseudo in inp.pseudos:
    print(pseudo)
```

Use  `set_vars` to set the value of several variables with a single call:

```{code-cell} ipython3
inp.set_vars(ecut=8, paral_kgb=0)
```

`AbinitInput` is a dict-like object, hence one can test for the presence of a variable in the input:

```{code-cell} ipython3
"ecut" in inp
```

To list all the variables that have been defined, use:

```{code-cell} ipython3
list(inp.keys())
```

To access the value of a particular variable use the syntax:

```{code-cell} ipython3
inp["ecut"]
```

To iterate over keywords and values:

```{code-cell} ipython3
for varname, varvalue in inp.items():
    print(varname, "-->", varvalue)
```

Use lists, tuples or numpy arrays when Abinit expects arrays

```{code-cell} ipython3
inp.set_vars(kptopt=1,
             ngkpt=[2, 2, 2],
             nshiftk=2,
             shiftk=[0.0, 0.0, 0.0, 0.5, 0.5, 0.5]  # 2 shifts in one list
            )

# It is possible to use strings but use them only for special cases such as:
inp["istwfk"] = "*1"
inp
```

If you mistype the name of the variable, `AbinitInput` raises an error:

```{code-cell} ipython3
try:
    inp.set_vars(perl=0)
except Exception as exc:
    print(exc)
```

<div class="alert alert-warning" role="alert">
The AbinitInput is a mutable object so changing it will aftect all the references
to the object.
See http://docs.python-guide.org/en/latest/writing/gotchas/ for further info
</div>

```{code-cell} ipython3
a = {"foo": "bar"}
b = a
c = a.copy()
a["hello"] = "world"
print("a dict:", a)
print("b dict:", b)
print("c dict:", c)
```
+++

The `set_structure` method sets the value of the ABINIT variables:

   * acell
   * rprim
   * ntypat
   * natom
   * typat
   * znucl
   * xred

It is always a good idea to set the structure immediately after the creation of `AbinitInput`
because several methods use this information to facilitate the specification of other variables.
For example, the `set_kpath` method uses the structure to generate the high-symmetry $k$-path for band structure calculations.

<div class="alert alert-warning">
`typat` must be consistent with the list of pseudopotentials passed to `AbinitInput`
</div>

+++

### Structure from dictionary:

```{code-cell} ipython3
structure = dict(
    ntypat=1,
    natom=2,
    typat=[1, 1],
    znucl=14,
    acell=3*[10.217],
    rprim=[[0.0,  0.5,  0.5],
           [0.5,  0.0,  0.5],
           [0.5,  0.5,  0.0]],
    xred=[[0.0 , 0.0 , 0.0],
          [0.25, 0.25, 0.25]]
)

inp = AbinitInput(structure, pseudos=abidata.pseudos("14si.pspnc"))
```

### Structure from file

+++

From a CIF file:

```{code-cell} ipython3
inp.set_structure(abidata.cif_file("si.cif"))
```

From a Netcdf file produced by ABINIT:

```{code-cell} ipython3
inp.set_structure(abidata.ref_file("si_scf_GSR.nc"))
```

Supported formats include:

   * *CIF*
   * *POSCAR/CONTCAR*
   * *CHGCAR*
   * *LOCPOT*
   * *vasprun.xml*
   * *CSSR*
   * *ABINIT netcdf files*
   * *pymatgen's JSON serialized structures*

+++

### From the Materials Project database:

```{code-cell} ipython3
# https://www.materialsproject.org/materials/mp-149/
inp.set_structure(abilab.Structure.from_mpid("mp-149"))
```

Remember to set the `PMG_MAPI_KEY` in ~/.pmgrc.yaml as described
[here](http://pymatgen.org/usage.html#setting-the-pmg-mapi-key-in-the-config-file).

+++

Note that you can avoid the call to `set_structure` if the `structure` argument is passed to `AbiniInput`:

```{code-cell} ipython3
AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))
```

+++ {"toc-hr-collapsed": true}

## Brillouin zone sampling

There are two different types of sampling of the BZ: homogeneous and high-symmetry k-path.
The later is mainly used for band structure calculations and requires the specification of:

   * kptopt
   * kptbounds
   * ndivsm

whereas the homogeneous sampling is needed for all the calculations in which
we have to compute integrals in the Brillouin zone e.g. total energy calculations, DOS, etc.
The $k$-mesh is usually specified via:

   * ngkpt
   * nshiftk
   * shiftk

+++

### Explicit $k$-mesh

```{code-cell} ipython3
inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))

# Set ngkpt, shiftk explicitly
inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=[0.0, 0.0, 0.0, 0.5, 0.5, 0.5])
```

### Automatic $k$-mesh

```{code-cell} ipython3
# Define a homogeneous k-mesh.
# nksmall is the number of divisions to be used to sample the smallest lattice vector,
# shiftk is automatically selected from an internal database.

inp.set_autokmesh(nksmall=4)
```

### High-symmetry $k$-path

```{code-cell} ipython3
# Generate a high-symmetry k-path (taken from an internal database)
# Ten points are used to sample the smallest segment,
# the other segments are sampled so that proportions are preserved.
# A warning is issued by pymatgen about the structure not being standard.
# Be aware that this might possibly affect the automatic labelling of the boundary k-points on the k-path.
# So, check carefully the k-point labels on the figures that are produced in such case.

inp.set_kpath(ndivsm=10)
```

## Utilities

Once the structure has been defined, one can compute the number of valence electrons with:

```{code-cell} ipython3
print("The number of valence electrons is: ", inp.num_valence_electrons)
```

If we need to change a particular (scalar) variable to generate inputs for convergence studies:

```{code-cell} ipython3
# When using a non-integer step, such as 0.1, the results will often not
# be consistent.  It is better to use ``linspace`` for these cases.
# See also numpy.arange and numpy.linspace

ecut_inps = inp.arange("ecut", start=2, stop=5, step=2)

print([i["ecut"] for i in ecut_inps])
```

```{code-cell} ipython3
tsmear_inps = inp.linspace("tsmear", start=0.001, stop=0.003, num=3)
print([i["tsmear"] for i in tsmear_inps])
```

## Invoking Abinit with AbinitInput

Once you have an `AbinitInput`, you can call Abinit to get useful information
or simply to validate the input file before running the calculation.
All the method that invoke Abinit starts with the `abi` prefix
followed by a verb e.g. `abiget` or `abivalidate`.

```{code-cell} ipython3
inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))

inp.set_vars(ecut=-2)
inp.set_autokmesh(nksmall=4)

v = inp.abivalidate()
if v.retcode != 0:
    # If there is a mistake in the input, one can acces the log file of the run with the log_file object
    print("".join(v.log_file.readlines()[-10:]))
```

Let's fix the problem with the negative ecut and rerun abivalidate!

```{code-cell} ipython3
inp["ecut"] = 2
inp["toldfe"] = 1e-10
v = inp.abivalidate()
if v.retcode == 0:
    print("All ok")
else:
    print(v)
```

At this point, we have a valid input file and we can get the k-points in the irreducible zone with:

```{code-cell} ipython3
ibz = inp.abiget_ibz()
print("number of k-points:", len(ibz.points))
print("k-points:", ibz.points)
print("weights:", ibz.weights)
print("weights are normalized to:", ibz.weights.sum())
```

We can also call the Abinit spacegroup finder with:

```{code-cell} ipython3
abistruct = inp.abiget_spacegroup()
print("spacegroup found by Abinit:", abistruct.abi_spacegroup)
```

To get the list of possible parallel configurations for this input up to 5 max_ncpus

```{code-cell} ipython3
inp["paral_kgb"] = 1
pconfs = inp.abiget_autoparal_pconfs(max_ncpus=5)
```

```{code-cell} ipython3
print("best efficiency:\n", pconfs.sort_by_efficiency()[0])
print("best speedup:\n", pconfs.sort_by_speedup()[0])
```

To get the list of irreducible phonon perturbations at Gamma (Abinit notation)

```{code-cell} ipython3
inp.abiget_irred_phperts(qpt=(0, 0, 0))
```

## Multiple datasets

Multiple datasets are handy when you have to generate several input files sharing several common
variables e.g. the crystalline structure, the value of ecut etc...
In this case, one can use the `MultiDataset` object that is essentially
a list of `AbinitInput` objects. Note however that `Abipy` workflows do not support input files with more than one dataset.

```{code-cell} ipython3
# A MultiDataset object with two datasets (a.k.a. AbinitInput)
multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                            pseudos="14si.pspnc", pseudo_dir=abidata.pseudo_dir, ndtset=2)

# A MultiDataset is essentially a list of AbinitInput objects
# with handy methods to perform global modifications.
# i.e. changes that will affect all the inputs in the MultiDataset
# For example:
multi.set_vars(ecut=4)

# is equivalent to
#
#   for inp in multi: inp.set_vars(ecut=4)
#
# and indeed:

for inp in multi:
    print(inp["ecut"])
```

```{code-cell} ipython3
# To change the values in a particular dataset use:
multi[0].set_vars(ngkpt=[2,2,2], tsmear=0.004)
multi[1].set_vars(ngkpt=[4,4,4], tsmear=0.008)
```

To build a table with the values of "ngkpt" and "tsmear":

```{code-cell} ipython3
multi.get_vars_dataframe("ngkpt", "tsmear")
```

```{code-cell} ipython3
multi
```

<div class="alert alert-warning" role="alert">
Remember that in python we start to count from zero hence the first dataset has index 0.
</div>

+++

Calling set_structure on `MultiDataset` will set the structure of the inputs:

```{code-cell} ipython3
multi.set_structure(abidata.cif_file("si.cif"))

# The structure attribute of a MultiDataset returns a list of structures
# equivalent to [inp.structure for inp in multi]
print(multi.structure)
```

The function `split_datasets` return the list of `AbinitInput` stored in MultiDataset

```{code-cell} ipython3
inp0, inp1 = multi.split_datasets()
inp0
```

<div class="alert alert">
You can use `MultiDataset` to build your input files but remember that
`Abipy` workflows will never support input files with more than one dataset.
As a consequence, you should always pass an `AbinitInput` to the
AbiPy functions that are building `Tasks`, `Works` or `Flows`.
</div>

```{code-cell} ipython3
print("Number of datasets:", multi.ndtset)
```

```{code-cell} ipython3
# To create and append a new dataset (initialized from dataset number 1)
multi.addnew_from(1)
multi[-1].set_vars(ecut=42)
print("Now multi has", multi.ndtset, "datasets and the ecut in the last dataset is:",
      multi[-1]["ecut"])
```
