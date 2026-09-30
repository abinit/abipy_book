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

# AbinitInput object

Creating Abinit input files is one of the most repetitive and error-prone operations
one has to perform before running calculations.
To simplify this task, AbiPy provides the `AbinitInput` object,
a dict-like object that stores the Abinit variables and provides easy-to-use methods to automate the definition of other variables.

This notebook discusses how to create an `AbinitInput` and how to define the parameters of the calculation.
In the last part, we introduce the `MultiDataset` object, which is mainly designed to generate
multiple inputs sharing the same structure and the same list of pseudopotentials.
In another [notebook](./input_factories), we briefly discuss how to use factory functions
to automatically generate input objects for typical calculations.

See also the {{ AbinitInput }}

```{include} snippets/plotly_matplotlib_note.md
```

```{include} snippets/manager_note.md
```

## Creating an AbinitInput object

```{code-cell}
import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

import os
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.abilab import AbinitInput

# This line tells AbiPy we are running inside a notebook
abilab.enable_notebook()

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

To create an Abinit input, we must specify the paths of the pseudopotential files.
In this case, we have a pseudo named `14si.pspnc` located in the `abidata.pseudo_dir` directory.

```{code-cell}
inp = AbinitInput(structure=abidata.cif_file("si.cif"),
                  pseudos="14si.pspnc", pseudo_dir=abidata.pseudo_dir)
```

`print(inp)` prints our input.
At this point, the input is almost empty since only the structure and the pseudos have been specified.

```{code-cell}
print(inp)
```

Inside a Jupyter notebook, the input can be displayed in HTML format
with links to the official ABINIT documentation:

```{code-cell}
inp
```

The input *has* a structure object:

```{code-cell}
print(inp.structure)
```

and a list of `Pseudo` objects:

```{code-cell}
for pseudo in inp.pseudos:
    print(pseudo)
```

built by parsing the pseudopotential files passed to `AbinitInput`.

Use `set_vars` to set several variables with a single call:

```{code-cell}
inp.set_vars(ecut=8, paral_kgb=0)
```

Since `AbinitInput` is a dict-like object, we can test whether a variable is present in the input:

```{code-cell}
"ecut" in inp
```

To list all the variables that have been defined, use:

```{code-cell}
list(inp.keys())
```

To access the value of a particular variable, use:

```{code-cell}
inp["ecut"]
```

To iterate over names and values:

```{code-cell}
for varname, varvalue in inp.items():
    print(varname, "-->", varvalue)
```

Use lists, tuples or numpy arrays when Abinit expects arrays:

```{code-cell}
inp.set_vars(
    kptopt=1,
    ngkpt=[2, 2, 2],
    nshiftk=2,
    shiftk=[0.0, 0.0, 0.0, 0.5, 0.5, 0.5]  # 2 shifts in one list
)

# It is possible to use strings but use them only for special cases such as:

inp["istwfk"] = "*1"
inp
```

+++

If you mistype the name of a variable, `AbinitInput` raises an exception:

```{code-cell}
try:
    inp.set_vars(perl=0)
except Exception as exc:
    print(exc)
```

+++

```{warning}
`AbinitInput` is a mutable object, so changing it will affect all the references to the object.
See [this page](http://docs.python-guide.org/en/latest/writing/gotchas/) for further info.
```

```{code-cell}
a = {"foo": "bar"}
b = a
c = a.copy()
a["hello"] = "world"

print("a dict:", a)
print("b dict:", b)
print("c dict:", c)
```

+++

The `set_structure` method sets the values of the ABINIT variables:

* {{ acell }}
* {{ rprim }}
* {{ ntypat }}
* {{ natom }}
* {{ typat }}
* {{ znucl }}
* {{ xred }}

It is always a good idea to set the structure immediately after creating an `AbinitInput`,
because several methods use this information to facilitate the specification of other variables.
For instance, the `set_kpath` method uses the structure to generate the high-symmetry $k$-path
for band structure calculations.

```{warning}
{{ typat }} must be consistent with the list of pseudopotentials passed to `AbinitInput`.
```

+++

### Creating a structure from Abinit variables

A structure can be created in different ways.

The most explicit (and verbose) approach consists in passing a dictionary with ABINIT variables,
using Python lists (or lists of lists) when ABINIT expects a 1D (or a multidimensional) array:

```{code-cell}
si_struct = dict(
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

print(si_struct)
```

+++

If you already have a string with the Abinit variables, you can use the `from_abistring` class method:

```{code-cell}
lif_struct = abilab.Structure.from_abistring("""
acell      7.7030079150    7.7030079150    7.7030079150 Angstrom
rprim      0.0000000000    0.5000000000    0.5000000000
           0.5000000000    0.0000000000    0.5000000000
           0.5000000000    0.5000000000    0.0000000000
natom      2
ntypat     2
typat      1 2
znucl      3 9
xred       0.0000000000    0.0000000000    0.0000000000
           0.5000000000    0.5000000000    0.5000000000
""")

print(lif_struct)
```

This approach requires less typing, yet we still need to specify {{ntypat}}, {{znucl}} and {{typat}}.
Fortunately, *from_abistring* supports another Abinit-specific format in which the
fractional coordinates and the element symbols are specified via the *xred_symbols* variable.
In this case, {{ntypat}}, {{znucl}} and {{typat}} do not need to be specified, as they are automatically
computed from *xred_symbols*:

```{code-cell}
mgb2_struct = abilab.Structure.from_abistring("""

# MgB2 lattice structure.
natom   3
acell   2*3.086  3.523 Angstrom
rprim   0.866025403784439  0.5  0.0
       -0.866025403784439  0.5  0.0
        0.0                0.0  1.0

# Atomic positions
xred_symbols
 0.0  0.0  0.0 Mg
 1/3  2/3  0.5 B
 2/3  1/3  0.5 B
""")

print(mgb2_struct)
```

### Structure from file

From a CIF file:

```{code-cell}
inp.set_structure(abidata.cif_file("si.cif"))
```

From one of the netcdf files produced by ABINIT:

```{code-cell}
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

### From the Materials Project database

```{code-cell}
# https://www.materialsproject.org/materials/mp-149/
inp.set_structure(abilab.Structure.from_mpid("mp-149"))
```

Remember to set `PMG_MAPI_KEY` in ~/.pmgrc.yaml as described
[here](https://pymatgen.org/usage.html#setting-the-pmg_mapi_key-in-the-config-file).

+++

Note that you can skip the call to `set_structure` if the `structure` argument is passed to `AbinitInput`:

```{code-cell}
AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))
```

+++ {"toc-hr-collapsed": true}

## Brillouin zone sampling

There are two types of BZ sampling: homogeneous meshes and high-symmetry $k$-paths.
The latter is mainly used for band structure calculations and requires the specification of:

* {{kptopt}}
* {{kptbounds}}
* {{ndivsm}}

whereas the homogeneous sampling is needed for all calculations that
require integrals over the Brillouin zone, e.g. total energy, DOS, etc.
The $k$-mesh is usually specified via:

* {{ngkpt}}
* {{nshiftk}}
* {{shiftk}}

+++

### Explicit $k$-mesh

```{code-cell}
inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))

# Set ngkpt, shiftk explicitly
inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=[0.0, 0.0, 0.0, 0.5, 0.5, 0.5])
```

### Automatic $k$-mesh

```{code-cell}
# Define a homogeneous k-mesh.
# nksmall is the number of divisions used to sample the smallest lattice vector,
# shiftk is automatically selected from an internal database.

inp.set_autokmesh(nksmall=4)
```

### High-symmetry $k$-path

```{code-cell}
# Generate a high-symmetry k-path (taken from an internal database)
# Ten points are used to sample the smallest segment,
# the other segments are sampled so that proportions are preserved.
# A warning is issued by pymatgen about the structure not being standard.
# Be aware that this might possibly affect the automatic labelling of the boundary k-points on the k-path.
# So, check carefully the k-point labels on the figures that are produced in such case.

inp.set_kpath(ndivsm=10)
```

## Utilities

Once the structure has been defined, the number of valence electrons can be computed with:

```{code-cell}
print("The number of valence electrons is: ", inp.num_valence_electrons)
```

To generate inputs for convergence studies in which a particular (scalar) variable is changed, use:

```{code-cell}
# When using a non-integer step, such as 0.1, the results will often not
# be consistent.  It is better to use ``linspace`` for these cases.
# See also numpy.arange and numpy.linspace

ecut_inps = inp.arange("ecut", start=2, stop=5, step=2)

print([i["ecut"] for i in ecut_inps])
```

```{code-cell}
tsmear_inps = inp.linspace("tsmear", start=0.001, stop=0.003, num=3)
print([i["tsmear"] for i in tsmear_inps])
```

## Invoking Abinit with AbinitInput

Once you have an `AbinitInput`, you can call Abinit to get useful information
or simply to validate the input file before running the calculation.
All the methods that invoke Abinit start with the `abi` prefix
followed by a verb, e.g. `abiget` or `abivalidate`.

```{code-cell}
inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))

inp.set_vars(ecut=-2)
inp.set_autokmesh(nksmall=4)

v = inp.abivalidate()
if v.retcode != 0:
    # If there is a mistake in the input, one can access the log file of the run with the log_file object
    print("".join(v.log_file.readlines()[-10:]))
```

Let's fix the negative {{ecut}} and rerun `abivalidate`:

```{code-cell}
inp["ecut"] = 2
inp["toldfe"] = 1e-10

v = inp.abivalidate()

if v.retcode == 0:
    print("All ok")
else:
    print(v)
```

Now that we have a valid input file, we can get the $k$-points in the irreducible zone with:

```{code-cell}
ibz = inp.abiget_ibz()
print("number of k-points:", len(ibz.points))
print("k-points:", ibz.points)
print("weights:", ibz.weights)
print("weights are normalized to:", ibz.weights.sum())
```

We can also call the Abinit spacegroup finder with:

```{code-cell}
abistruct = inp.abiget_spacegroup()
print("spacegroup found by Abinit:", abistruct.abi_spacegroup)
```

To get the list of possible parallel configurations for this input with up to 5 CPUs ({{max_ncpus}}):

```{code-cell}
inp["paral_kgb"] = 1
pconfs = inp.abiget_autoparal_pconfs(max_ncpus=5)
```

```{code-cell}
print("best efficiency:\n", pconfs.sort_by_efficiency()[0])
print("best speedup:\n", pconfs.sort_by_speedup()[0])
```

To get the list of irreducible phonon perturbations at $\Gamma$ (Abinit notation):

```{code-cell}
inp.abiget_irred_phperts(qpt=(0, 0, 0))
```

## Multiple datasets

Multiple datasets are handy when you have to generate several input files sharing common
variables, e.g. the crystalline structure, the value of {{ecut}}, etc.
In this case, one can use the `MultiDataset` object, which is essentially
a list of `AbinitInput` objects.
Note, however, that AbiPy workflows do not support input files with more than one dataset.

```{code-cell}
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

```{code-cell}
# To change the values in a particular dataset use:
multi[0].set_vars(ngkpt=[2, 2, 2], tsmear=0.004)
multi[1].set_vars(ngkpt=[4, 4, 4], tsmear=0.008)
```

To build a table with the values of {{ngkpt}} and {{tsmear}}:

```{code-cell}
multi.get_vars_dataframe("ngkpt", "tsmear")
```

```{code-cell}
multi
```

```{warning}
Remember that Python counts from zero, hence the first dataset has index 0.
```

+++

Calling *set_structure* on a `MultiDataset` sets the structure of all the inputs:

```{code-cell}
multi.set_structure(abidata.cif_file("si.cif"))

# The structure attribute of a MultiDataset returns a list of structures
# equivalent to [inp.structure for inp in multi]
print(multi.structure)
```

The `split_datasets` method returns the list of `AbinitInput` objects stored in the `MultiDataset`:

```{code-cell}
inp0, inp1 = multi.split_datasets()
inp0
```

```{note}
You can use `MultiDataset` to build your input files, but remember that
AbiPy workflows will never support input files with more than one dataset.
As a consequence, you should always pass an `AbinitInput` to the
AbiPy functions that build `Tasks`, `Works` or `Flows`.
```

```{code-cell}
print("Number of datasets:", multi.ndtset)
```

```{code-cell}
# To create and append a new dataset (initialized from dataset number 1)
multi.addnew_from(1)
multi[-1].set_vars(ecut=42)
print("Now multi has", multi.ndtset, "datasets and the ecut in the last dataset is:",
      multi[-1]["ecut"])
```
