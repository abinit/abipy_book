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

Back to the main [Index](index.ipynb) <a id="top"></a>

+++

# Structure object

The `AbiPy` structure inherits from the `pymatgen` structure. 
One has therefore access to all the methods and tools already available in `pymatgen`.
In this notebook, we mainly focus on the extensions added by `AbiPy`. 
For the features provided by pymatgen, please consult the 
[official pymatgen documentation](http://pymatgen.org/usage.html#structures-and-molecules)

## Table of Contents
[[back to top](#top)]

- [Reading a structure from file](#Reading-a-structure-from-file)
- [Converting to other formats](#Converting-to-other-formats)
- [Getting information on the structure](#Getting-information-on-the-structure)
- [abistruct.py script](#abistruct.py)

```{code-cell} ipython3
from __future__ import division, print_function, unicode_literals

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
[[back to top](#top)]

+++

It is possible to initialize a crystalline structure from different file formats: 

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

```{code-cell} ipython3
# abidata.cif_file returns one of the CIF files shipped with AbiPy.
structure = Structure.from_file(abidata.cif_file("si.cif"))
print(structure)
```

To read the structure from a netcdf file:

```{code-cell} ipython3
structure = Structure.from_file(abidata.ref_file("si_nscf_GSR.nc"))

# Use to_string with verbose > 0 to get more info 
print(structure.to_string(verbose=1))
```

Use `to_abivars` to get the list of Abinit variables in a python dictionary:

```{code-cell} ipython3
structure.to_abivars()
```

and `abi_string` to get a string that can be used directly in the input file:

```{code-cell} ipython3
print(structure.abi_string)
```

To visualize the structure with matplotlib, use:

```{code-cell} ipython3
structure.plot();
```

The matplotlib version is minimalistic but it plays well with jupyter notebooks.
For a more advanced visualization we suggest using a specialized graphical applications.
Fortunately, one can invoke (already installed) external applications directly from AbiPy with e.g.

```{code-cell} ipython3
# structure.visualize("vesta")
```

To get a structure from the materials project database 
(https://www.materialsproject.org ), use:

```{code-cell} ipython3
# You can pass the api_key or set the env variable PMG_MAPI_KEY in your ~/.pmgrc.yaml files.
si2_mp = Structure.from_mpid("mp-149", api_key=None)
print(si2_mp)
```

In some cases, we have multiple structures and we need to compare the lattice parameters. 
Use `dataframes_from_structures` to build a pandas DataFrame:

```{code-cell} ipython3
dfs = abilab.dataframes_from_structures([structure, si2_mp], index=["CIF", "MP"])
```

then we can compare the lattice parameters with:

```{code-cell} ipython3
dfs.lattice
```

Note that all AbiPy robots have this feature built-in. 
Sometimes it is much easier to build a robot directly from files 
and then compare the structures with e.g. `robot.get_lattice_dataframe()`.

+++

## Converting to other formats
[[back to top](#top)]

+++

Use `structure.convert(format)` to get the string representation in the new format:

```{code-cell} ipython3
for fmt in ["cif", "POSCAR", "qe"]:
    print((" Abinit --> %s " % fmt).center(80, "*"))
    print(structure.convert(fmt=fmt))
```

## Getting information on the structure
[[back to top](#top)]

```{code-cell} ipython3
print(structure.reciprocal_lattice)
```

```{code-cell} ipython3
structure.reciprocal_lattice.matrix.T @ structure.lattice.matrix / (2 * np.pi)
```

```{code-cell} ipython3
# List of high-symmetry k-points.
print(structure.hsym_kpoints)
```

The method `calc_ksampling` allows one to get an efficient sampling of the Brillouin zone 
by just specifying the number of divisions to be used for the smallest lattice vector of the reciprocal lattice:

```{code-cell} ipython3
pprint(structure.calc_ksampling(nksmall=10))
```

To get the recommended high symmetry $k$-path in reduced coordinates:

```{code-cell} ipython3
structure.calc_kptbounds()
```

The high-symmetry q-path is automatically selected assuming
the structure fulfills the convention described in [Setyawan2010](https://doi.org/10.1016/j.commatsci.2010.05.010)

+++

To visualize the Brillouin zone with matplotlib:

```{code-cell} ipython3
structure.plot_bz();
```

To get the number of valence electrons for a given set of pseudopotentials: 

```{code-cell} ipython3
structure.num_valence_electrons(pseudos=abidata.pseudos("14si.pspnc"))
```

To visualize the X-ray diffraction plot with pymatgen XRDCalculator

```{code-cell} ipython3
structure.plot_xrd();
```

## abistruct.py 
[[back to top](#top)]

`abistruct.py` provides a handy command line interface to operate on structure objects 
constructed from external files. 
There are several options available as well an interface to the [materials project](http://materialsproject.org/)
and the [COD](http://www.crystallography.net/cod/) database.

```{code-cell} ipython3
!abistruct.py --help
```

Back to the main [Index](index.ipynb)
