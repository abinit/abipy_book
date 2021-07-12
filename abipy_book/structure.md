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

# Structure object

The `AbiPy` structure inherits from the `pymatgen` structure.
One has therefore access to all the methods and tools already available in `pymatgen`.
In this notebook, we mainly focus on the extensions added by `AbiPy`.
For the features provided by pymatgen, please consult the
[official pymatgen documentation](http://pymatgen.org/usage.html#structures-and-molecules).

{{ Structure }} and {{ pymatgen_Structure}}

```{code-cell}
import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

# Import abipy modules
from abipy import abilab
from abipy.abilab import Structure
import abipy.data as abidata

# Useful tools we'll need later on.
from pprint import pprint
import numpy as np

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

## Reading a structure from file

It is possible to initialize a structure object from different file formats:

   * CIF
   * POSCAR/CONTCAR
   * CHGCAR
   * LOCPOT,
   * vasprun.xml
   * CSSR
   * ABINIT Netcdf files
   * pymatgen's JSON serialized structures

Note, in particular, that one can initialize the structure from the netcdf files
produced by Abinit (`GSR.nc`, `WFK.nc`, etc) as well as output files in text format
such as the Abinit input/output files or even the DDB file.

To initialize the structure from a CIF file use the `from_file` method:

```{code-cell}

structure = Structure.from_file(abidata.cif_file("si.cif"))
print(structure)
```

```{include} abidata_note.md
```

To read the structure from an Abinit netcdf file, use:

```{code-cell}
structure = Structure.from_file(abidata.ref_file("si_nscf_GSR.nc"))

print(structure.to_string(verbose=1))  # Use to_string with verbose > 0 to get more info
```

Use `to_abivars` to get a python dictionary with the list of Abinit variables.

```{code-cell}
structure.to_abivars()
```

and the `abi_string` property to get a string that can be used directly in the input file:

```{code-cell}
print(structure.abi_string)
```

To visualize the structure with matplotlib, use:

```{code-cell}
structure.plot();
```

The matplotlib version is minimalistic but it plays well with jupyter notebooks.
For a more advanced visualization we suggest using a specialized graphical applications.
Fortunately, one can invoke external applications directly from AbiPy with e.g.

```{code-cell}
# structure.visualize("vesta")
```

provided VESTA is already installed on your machine and the binary can be found in  **$PATH**.

To get a structure from the [materials project database](https://www.materialsproject.org), use:

```{code-cell}
# You can pass the api_key or set the env variable PMG_MAPI_KEY in your ~/.pmgrc.yaml files.
si2_mp = Structure.from_mpid("mp-149", api_key=None)
print(si2_mp)
```

In some cases, we have multiple structures and we need to compare the lattice parameters.
Use `dataframes_from_structures` to build a pandas DataFrame:

```{code-cell}
dfs = abilab.dataframes_from_structures([structure, si2_mp], index=["CIF", "MP"])
```

then we can compare the lattice parameters with:

```{code-cell}
dfs.lattice
```

Note that all AbiPy robots have this feature built-in.
Sometimes it is much easier to build a robot directly from files
and then compare the structures with e.g. `robot.get_lattice_dataframe()`.


## Converting to other formats

Use `structure.convert(format)` to get the string representation in the new format:

```{code-cell}
for fmt in ("cif", "POSCAR", "qe"):
    print((" Abinit --> %s " % fmt).center(80, "*"))
    print(structure.convert(fmt=fmt))
```

## Getting info on the structure

```{code-cell}
print(structure.reciprocal_lattice)
```

```{code-cell}
structure.reciprocal_lattice.matrix.T @ structure.lattice.matrix / (2 * np.pi)
```

```{code-cell}
# List of high-symmetry k-points.
print(structure.hsym_kpoints)
```

The method `calc_ksampling` allows one to get an efficient sampling of the Brillouin zone
by just specifying the number of divisions to be used for the smallest lattice vector of the reciprocal lattice:

```{code-cell}
pprint(structure.calc_ksampling(nksmall=10))
```

To get the recommended high symmetry $k$-path in reduced coordinates:

```{code-cell}
structure.calc_kptbounds()
```

The high-symmetry **q**-path is automatically selected assuming
the structure fulfills the convention described in [Setyawan2010](https://doi.org/10.1016/j.commatsci.2010.05.010)

+++

To visualize the Brillouin zone with matplotlib, use:

```{code-cell}
structure.plot_bz();
```

For the plotly version, use:

```{code-cell}
structure.plotly_bz();
```

```{note}
The name of the plotly method (if implemented) is obtained by replacing the `plot` verb with `plotly`.
```

To get the number of valence electrons for a given set of pseudopotentials:

```{code-cell}
structure.num_valence_electrons(pseudos=abidata.pseudos("14si.pspnc"))
```

To visualize the X-ray diffraction plot with pymatgen XRDCalculator, use:

```{code-cell}
structure.plot_xrd();
```

+++

## The `abistruct.py` script

The {{ abistruct }} script provides a handy command line
interface to operate on structure objects constructed from external files.
There are several options available as well an interface to the {{ materials project }}
and the {{ COD }} database.

To obtain the list of available commands, use:

```{code-cell}
!abistruct.py --help
```

## Creating a GUI inside a notebook

Several AbiPy objects provide a `get_panel` method that allows one to create a {{ panel }} GUI
exposing some of the underlying AbiPy methods.
Similar capabilities are also available via the {{ abipygui }} web app.

To build a panel GUI for a given structure use:

```{code-cell}
abilab.abipanel()
structure.get_panel()
```
