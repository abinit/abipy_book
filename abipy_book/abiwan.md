---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python [default]
  language: python
  name: python3
---

# The ABIWAN file (wannier90)

This notebook shows how to use AbiPy to analyze the output files
produced by [wannier90](http://www.wannier.org/) and how to use the `ABIWAN.nc` file
produced by Abinit to interpolate band energies.

As usual, one can use:

    abiopen.py FILE_ABIWAN.nc

with the `--expose` or the `--print` option for a command line interface
and `--notebook` to generate a jupyter notebook.

Note: The code in this notebook requires abinit >= 8.9 and abipy >= 0.6


## Table of Contents

- [How to analyze the WOUT file](#How-to-analyze-the-WOUT-file)
- [Using ABIWAN.nc to interpolate band energies](#Using-ABIWAN.nc-to-interpolate-band-energies)

Let's start by importing the basic modules needed for this tutorial.

```{code-cell} ipython3
import os

import warnings
warnings.filterwarnings("ignore")  # Ignore warnings

from abipy import abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook
import abipy.data as abidata

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

## How to analyze the WOUT file

+++

Use `abiopen` to open a `wout` file (the main output file produced by wannier90):

```{code-cell} ipython3
filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "example01_gaas.wout")

wout = abilab.abiopen(filepath)
print(wout)
```

To plot the convergence of the wannier cycle use:

```{code-cell} ipython3
wout.plot();
```

To plot the evolution of the Wannier centers and spread, use:

```{code-cell} ipython3
wout.plot_centers_spread();
```

<div class="alert alert-info" role="alert">
Alternatively one can use `abiopen.py FILE_MDF.nc -nb` to generate a jupyter notebook directly from the terminal
or `abiopen.py FILE_MDF.nc -e -sns` to produce matplotlib plots automatically.
</div>

+++

## Using ABIWAN.nc to interpolate band energies


`ABIWAN.nc` is a netcdf file produced by Abinit after having called *wannier90* in library mode.
The file contains the unitary transformation and other important parameters associated to the calculations.
This file can be read by AbiPy and can be used to interpolate band energies with the wannier method.

As usual, use `abiopen` to open the file:

```{code-cell} ipython3
filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "tutoplugs_tw90_4", "tw90_4o_DS3_ABIWAN.nc")
abiwan = abilab.abiopen(filepath)
print(abiwan)
```

To plot the matrix elements of the KS Hamiltonian in real space in the Wannier Gauge, use:

```{code-cell} ipython3
abiwan.hwan.plot(title="Matrix elements in real space");
```

To interpolate the KS energies along a high-symmetry k-path and construct
a new `ElectronBands` object, use:

```{code-cell} ipython3
ebands_kpath = abiwan.interpolate_ebands()
```

```{code-cell} ipython3
ebands_kpath.plot(title="Wannier-interpolated");
```

If you need an IBZ sampling instead of a k-path, for instance a 36x36x36 k-mesh, use:

```{code-cell} ipython3
ebands_kmesh = abiwan.interpolate_ebands(ngkpt=(36, 36, 36))
```

As we are dealing with AbiPy objects, we can easily reuse the AbiPy API to plot bands + DOS:

```{code-cell} ipython3
ebands_kpath.plot_with_edos(ebands_kmesh.get_edos(), title="Wannier-interpolated bands and DOS");
```

We can also compare an ab-initio band structure with the Wannier-interpolated results.
This is useful to understand if our wannier functions are well localized and if the
k-mesh used with wannier90 is dense enough.

In this case, it is just a matter of passing the path to the netcdf file
containing the ab-initio band structure to the `get_plotter_from_ebands` method of `abiwan`.
The function interpolates the band energies using the k-path found in the netcdf file
and returns a plotter object:

```{code-cell} ipython3
import abipy.data as abidata
gsr_path = abidata.ref_file("si_nscf_GSR.nc")

plotter = abiwan.get_plotter_from_ebands(gsr_path)
```

Then we call `combiplot` to plot the two band structures on the same figure:

```{code-cell} ipython3
plotter.combiplot();
```

As we can see, the interpolated band structures is not completely on top of the ab-initio
results. To improve the agreement we should try to reduced the spread and/or increase
the density of the k-mesh used in the wannierization procedure.

+++
