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

# The HIST.nc file (relaxation/MD)

The `HIST.nc` file contains the history of structural relaxations or molecular dynamics calculations.
One can use the `abiopen` function provide by `abilab` to open the file and generate an instance of `HistFile`.
Alteratively, one can use the `abiopen.py` script to open the file inside the shell with the syntax:

    abiopen.py out_HIST.nc

This command will start the ipython interpreter so that one can interact directly
with the `HistFile` object (named `abifile` inside ipython).

To generate a jupyter notebook use:

    abiopen.py out_HIST.nc -nb

For a quick visualization of the data, use the `--expose` option:

    abiopen.py out_HIST.nc -e

```{code-cell}
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

```{code-cell}
hist = abilab.abiopen(abidata.ref_file("sic_relax_HIST.nc"))
print("Number of iterations performed:", hist.num_steps)
```

`hist.structures` is the list of structure objects at the different iteration steps.
`hist.etotals` is a numpy array with the total energies in eV associated to the different steps.

```{code-cell}
for struct, etot in zip(hist.structures, hist.etotals):
    print("Volume:", struct.volume,", Etotal:", etot)
```

To get the last structure stored in the `HIST.nc` file:

```{code-cell}
print(hist.final_structure)
```

To plot the evolution of the structural parameters with `matplotlib`:

```{code-cell}
hist.plot(tight_layout=True);
```

```{code-cell}
hist.plotly();
```

To plot the total energies at the different iterations steps:

```{code-cell}
hist.plot_energies();
```

```{code-cell}
hist.plotly_energies();
```

## Converting to other formats

Use `to_xdatcar` to get a XDATCAR pymatgen object (useful to interface AbiPy with other pymatgen tools)

```{code-cell}
# hist.write_xdatcar writes a XDATCAR file
xdatcar = hist.to_xdatcar()
print(xdatcar)
```

