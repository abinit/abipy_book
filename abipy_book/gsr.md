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

# The GSR file (Ground-State)

In this notebook we discuss how to plot the electron band structures and the density of states (DOS)
using the GSR netcdf files produced by Abinit.

For the tutorial, we will use the netcdf files shipped with AbiPy.
The function `abidata.ref_file` returns the absolute path of the reference file.
In your scripts, you have to replace `data.ref_file("abipy_filename")` with a string giving
the location of your netcdf file.

Alteratively, one can use the `abiopen.py` script to open the file inside the shell with the syntax:

    abiopen.py out_GSR.nc

This command will start the ipython interpreter so that one can interact directly
with the `GsrFile` object (named `abifile` inside ipython).
To generate a jupyter notebook use:

    abiopen.py out_GSR.nc -nb

For a quick visualization of the data, use the `--expose` options:

    abiopen.py out_GSR.nc -e

```{include} snippets/abiopen_note.md
```

```{code-cell} 
import warnings
warnings.filterwarnings("ignore")  # Ignore warnings

from abipy import abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook
# Import abipy reference data.
import abipy.data as abidata

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

## The GSR File

The `GSR` file (mnemonics: Ground-State Results) is a netcdf file with the
results produced by SCF or NSCF ground-state calculations
(band energies, forces, energies, stress tensor).
To open a `GSR` file, use the `abiopen` function defined in `abilab`:

```{code-cell} 
gsr = abilab.abiopen(abidata.ref_file("si_scf_GSR.nc"))
```

The gsr object has a `Structure`:

```{code-cell} 
print(gsr.structure)
```

and an `ElectronBands` object with the band energies, the occupation factors, the list of k-points:

```{code-cell} 
print(gsr.ebands)
```

```{important}
In python we start to count from zero, thus the first band has index 0 and the first spin is 0
AbiPy uses the same convention so be very careful when specifying band, spin or k-point indices.
```

A GSR file produced by a **self-consistent run**, contains the values of the total energy, the forces,
and the stress tensor at the end of the SCF cycle:

```{code-cell} 
print("energy:", gsr.energy, "pressure:", gsr.pressure)
```

To get a summary of the most important results:

```{code-cell} 
print(gsr)
```

The different contributions to the total energy are stored in a dictionary:

```{code-cell} 
print(gsr.energy_terms)
```

At this point, we don't need this file anymore so we close it with:

```{code-cell} 
gsr.close()
```

```{warning}
The gsr maintains a reference to the underlying netcdf file hence one should
call `gsr.close()` to release the resource when we don't need it anymore.
Python will do it automatically if you use `abiopen` and the `with` context manager.

Note that we don't always follow this rule inside the jupyter notebook to maintain the
code readable but you should definitively close all your files, especially when
writing code that may be running for hours or even more.
```

## Plotting band structures

Let's open the GSR file produced by a NSCF calculation done on a high-symmetry k-path
and extract the electronic band structure.

A warning is issued by pymatgen about the structure not being standard.
Be aware that this might possibly affect the automatic labelling of the boundary k-points on the k-path.
So, check carefully the k-point labels on the figures that are produced in such case.
In the present case, the labelling is correct.

```{code-cell} 
with abilab.abiopen(abidata.ref_file("si_nscf_GSR.nc")) as nscf_gsr:
    ebands_kpath = nscf_gsr.ebands
```

Now we can plot the band energies with *matplotlib*:

```{code-cell} 
# The labels for the k-points are found in an internal database.
ebands_kpath.plot(with_gaps=True, title="Silicon band structure");
```

Alternatively, one can use the optional argument `klabels` to define the mapping
`reduced_coordinates --> name of the k-point` and pass it to the plot method

```{code-cell} 
klabels = {
    (0.5, 0.0, 0.0): "L",
    (0.0, 0.0, 0.0): "$\Gamma$",
    (0.0, 0.5, 0.5): "X"
}

# ebands_kpath.plot(title="User-defined k-labels", band_range=(0, 5), klabels=klabels);
```

For the plotly version, use:

```{code-cell} 
ebands_kpath.plotly(with_gaps=True, title="Silicon band structure with plotly");
```

```{code-cell} 
abilab.abipanel()
gsr.get_panel()
```

Let's have a look at our k-points by calling `kpoints.plot()`

```{code-cell} 
ebands_kpath.kpoints.plotly();
```

and the crystalline structure with:

```{code-cell} 
ebands_kpath.structure.plot();
```

```{note}
The same piece of code works if you replace the `GSR.nc` file with e.g. a `WFK.nc` file in netcdf format
(actually any netcdf file with an ebands object).
The main advantage of the `GSR` file is that it is lightweight (no wavefunctions).
```

## DOS with the Gaussian technique

Let's use the eigenvalues and the k-point weights stored in `gs_ebands` to
compute the DOS with the Gaussian method.
The method is called without arguments so we use **default values**
for the *broadening* and the *step* of the linear mesh.

```{code-cell} 
 with abilab.abiopen(abidata.ref_file("si_scf_GSR.nc")) as scf_gsr:
    ebands_kmesh = scf_gsr.ebands

edos = ebands_kmesh.get_edos()
print(edos)
```

```{code-cell} 
edos.plotly();
```

```{code-cell} 
print("[ebands_kmesh] is_ibz:", ebands_kmesh.kpoints.is_ibz, "is_kpath:", ebands_kmesh.kpoints.is_path)
print("[ebands_kpath] is_ibz:", ebands_kpath.kpoints.is_ibz, "is_kpath:", ebands_kpath.kpoints.is_path)
```

```{warning}
The DOS requires a homogeneous $k$-sampling of the BZ. Abipy will raise an exception if you try
to compute the DOS with a k-path.
```

To plot bands and DOS on the same figure:

```{code-cell} 
ebands_kpath.plotly_with_edos(edos, with_gaps=True);
```

To plot the DOS and the integrated DOS (IDOS), use:

```{code-cell} 
edos.plotly_dos_idos();
```

The gaussian broadening can significantly change the overall shape of the DOS.
If accurate values are needed (e.g. the DOS at the Fermi level in metals),
one should perform an accurate convergence study with respect to the k-point mesh.
Here we show how compute the DOS with different values of the gaussian smearing
for fixed k-point sampling and plot the results:

```{code-cell} 
# Compute the DOS with the Gaussian method and different values of the broadening
widths = [0.1, 0.2, 0.3, 0.4]

edos_plotter = ebands_kmesh.compare_gauss_edos(widths, step=0.1)
```

To plot the results on the same figure, use:

```{code-cell} 
edos_plotter.combiplot(title="e-DOS as function of the Gaussian broadening");
```

while `gridplot` generates a grid of subplots:

```{code-cell} 
edos_plotter.gridplot();
```

## Joint density of states

This example shows how to plot the different contributions to the electronic joint density of states of Silicon.
Select the valence and conduction bands to be included in the JDOS. Here we include valence bands from 0 to 3 and the first conduction band (4).

```{code-cell} 
vrange = range(0,4)
crange = range(4,5)

# Plot data
ebands_kmesh.plot_ejdosvc(vrange, crange);
```

```{code-cell} 
ebands_kpath.plot_transitions(omega_ev=3.0);
```

```{code-cell} 
ebands_kmesh.plot_transitions(omega_ev=3.0);
```

## Plotting the Fermi surface

```{code-cell} 
with abilab.abiopen(abidata.ref_file("mgb2_kmesh181818_FATBANDS.nc")) as fbnc_kmesh:
    mgb2_ebands = fbnc_kmesh.ebands

    # Build ebands in full BZ.
    mgb2_eb3d = mgb2_ebands.get_ebands3d()
```

There are three bands crossing the Fermi level of $MgB_2$ (band 2, 3, 4):

```{code-cell} 
mgb2_ebands.boxplot();
```

Let's use matplotlib to plot the isosurfaces corresponding to the Fermi level (default):

```{code-cell} 
# Warning: requires skimage package, rendering could be slow.
mgb2_eb3d.plot_isosurfaces();
```

```{code-cell} 
#ebands_kpath.plot_scatter3d(band=3);
#ebands_kmesh.plot_scatter3d(band=3);
```

## Analyzing multiple GSR files with robots

TODO

```{note}
Robots can also be constructed from the command line with: abicomp.py gsr FILES
```

