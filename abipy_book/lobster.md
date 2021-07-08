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

# Lobster output files

This example shows how to analyze the output files
produced by [Lobster](http://schmeling.ac.rwth-aachen.de/cohp)

Use

    abiopen.py FILE

with the `--expose` or `--print` for a command line interface
and `--notebook` to generate a jupyter notebook from a lobster `FILE`.

Note: The code in this notebook requires abipy >= 0.6

Let's start by importing the basic modules needed for this tutorial.

```{code-cell} ipython3
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

## How to analyze the COHPCAR file

```{code-cell} ipython3
# Path to one of the reference file shipped with AbiPy
import os
dirpath = os.path.join(abidata.dirpath, "refs", "lobster_gaas")
filename = os.path.join(dirpath, "GaAs_COHPCAR.lobster.gz")

# Open the COHPCAR.lobster file (same API for COOPCAR.lobster)
cohp_file = abilab.abiopen(filename)
print(cohp_file)
```

To plot the COHP averaged over all atom pairs specified:

```{code-cell} ipython3
cohp_file.plot(title="GaAs COHP");
```

To plot the integrated COHP averaged over all atom pairs:

```{code-cell} ipython3
cohp_file.plot(what="i", title="GaAs integrated COHP");
```

To plot the total overlap for all sites listed in `from_site_index`

```{code-cell} ipython3
cohp_file.plot_site_pairs_total(from_site_index=[0, 1], title="COHP total overlap for site index 0");
```

To plot partial crystal orbital projections for all sites listed in `from_site_index`:

```{code-cell} ipython3
cohp_file.plot_site_pairs_partial(from_site_index=[0, 1],
                                  title="COHP with orbital projections from site index 0",
                                  fontsize=6, tight_layout=True);
```

```{code-cell} ipython3
#cohp_file.plot_average_pairs(with_site_index=[0]);
```

Use `abiopen` to open the MDF:

+++

## How to analyze the ICOHPLIST file

```{code-cell} ipython3
# Path to one of the AbiPy file
dirpath = os.path.join(abidata.dirpath, "refs", "lobster_gaas")
filename = os.path.join(dirpath, "GaAs_ICOHPLIST.lobster.gz")

# Open the ICOHPCAR.lobster file.
icohp_file = abilab.abiopen(filename)
print(icohp_file)
```

## How to analyze the DOSCAR file

```{code-cell} ipython3
dirpath = os.path.join(abidata.dirpath, "refs", "lobster_gaas")
filename = os.path.join(dirpath, "GaAs_DOSCAR.lobster.gz")

# Open the ICOHPCAR.lobster file.
doscar = abilab.abiopen(filename)
print(doscar)
```

```{code-cell} ipython3
doscar.plot();
```

```{code-cell} ipython3
doscar.plot_pdos_site(site_index=[0, 1]);
```

## Analyzing all Lobster output files with LobsterAnalyzer

Let's assume we have a directory with lobster output files
for COOP, COHP, DOS and we need to produce plots showing all these results altogether.
In this case, one can use the `LobsterAnalyzer` object and initialize it from the directory
containing the output files.

```{code-cell} ipython3
dirpath = os.path.join(abidata.dirpath, "refs", "lobster_gaas")

# Open the all the lobster files produced in directory dirpath
# with the (optional) prefix GaAs_
lobana = abilab.LobsterAnalyzer.from_dir(dirpath, prefix="GaAs_")
print(lobana.to_string(verbose=1))
```

To plot COOP + COHP + DOS, use:

```{code-cell} ipython3
lobana.plot(title="COOP + COHP + DOS");
```

To plot COHP for all sites in from_site_index and Lobster DOS:

```{code-cell} ipython3
lobana.plot_coxp_with_dos(from_site_index=[0, 1]);
```

```{code-cell} ipython3
# Plot orbital projections.
lobana.plot_coxp_with_dos(from_site_index=[0], with_orbitals=True);
```

```{code-cell} ipython3
#lobana.plot_with_ebands(ebands="out_GSR.nc")
```

<div class="alert alert-info" role="alert">
For a command line interface, use: `abiview.py lobster .`.
Use the `--expose` option to generate plots automatically.
</div>
