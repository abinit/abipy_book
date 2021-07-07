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

Back to the main [Index](index.ipynb) <a id="top"></a>

+++

# Postprocessing tools for Bethe-Salpeter calculations

+++

The Bethe-Salpeter code saves the optical spectra in the `MDF.nc` file.
This notebook explains how to use the AbiPy API to analyze the results.

## Table of Contents
[[back to top](#top)]

- [How to analyze a single MDF file](#How-to-analyze-a-single-MDF-file)
- [Analyzing multiple MDF files with robots](#Analyzing-multiple-MDF-files-with-robots)

Let's start by importing the basic modules we will need for this tutorial.

```{code-cell} ipython3
# Use this at the beginning of your script so that your code will be compatible with python3
from __future__ import print_function, division, unicode_literals

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

## How to analyze a single MDF file  
[[back to top](#top)]

+++

Use `abiopen` to open the MDF:

```{code-cell} ipython3
mdf_file = abilab.abiopen(abidata.ref_file("tbs_4o_DS2_MDF.nc"))
print(mdf_file)
```

To plot the (averaged) imaginary part of the macroscopic dielectric function (MDF)
between 2 and 5 eV use:

```{code-cell} ipython3
mdf_file.plot_mdfs(title="Si absorption spectrum: EXC vs RPA", xlims=(2, 5));
```

To select the MDF computed for the first q-point, use

```{code-cell} ipython3
mdf_file.plot_mdfs(title="Im(Mdf) at the first q-point", qpoint=0, xlims=(2, 5));
```

* EXC: MDF with excitonic effects included
* KS-RPA: MDF computed with KS eigenvalues
* GW-RPA: MDF computed at the RPA level with KS + scissors operator

+++

To plot the (averaged) real part of the MDF:

```{code-cell} ipython3
mdf_file.plot_mdfs(cplx_mode="re", title="Real part of MDF: EXC vs RPA", xlims=(2, 5));
```

<div class="alert alert-info" role="alert">
Alternatively one can use `abiopen.py FILE_MDF.nc -nb` to generate a jupyter notebook directly from the terminal
or `abiopen.py FILE_MDF.nc -e -sns` to produce matplotlib plots automatically.
</div>

+++

## Analyzing multiple MDF files with robots
[[back to top](#top)]

To analyze the converge of the optical spectra, we can use the MdfRobot.
Let's build our robot from a list of MDF.nc files:

```{code-cell} ipython3
paths = abidata.ref_files("si_444_MDF.nc", "si_666_MDF.nc", "si_888_MDF.nc")
robot = abilab.MdfRobot.from_files(paths)
print(robot)
```

```{code-cell} ipython3
plotter = robot.get_multimdf_plotter()
```

To analyze the convergence of the (averaged) MDFs:

```{code-cell} ipython3
plotter.plot();
```

It is also possible to analyze the converge of the MDF for the different q-directions with:

```{code-cell} ipython3
plotter.plot(qview="all");
```

<div class="alert alert-info" role="alert">
Robots can also be constructed from the command line with: `abicomp.py mdf FILES`.
Use the `--expose` option to generate plots automatically.
</div>

Back to the main [Index](index.ipynb)
