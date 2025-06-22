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

# The SIGRES file (GW)

This notebook explains how to use AbiPy and `matplotlib` to visualize the results produced by the GW code.
The self-energy code ({{optdriver}} 4) saves the final results in the `SIGRES.nc` file
while the screening code ({{optdriver}} 3) stores the inverse dielectric matrix in the `SCR.nc` file.

Let's start by importing the basic modules we will need for this tutorial.

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

## How to visualize QP corrections

As usual, we start by opening the netcdf file with abiopen:

```{code-cell} ipython3
sigres = abilab.abiopen(abidata.ref_file("tgw1_9o_DS4_SIGRES.nc"))
print(sigres)
```

Let's have a look at the KS energies used to compute the Green's function $G_0$, the RPA screening $W_0$
and the $G_0W_0$ self-energy:

```{code-cell} ipython3
sigres.ebands.plot();
```

The SIGRES file contains the KS as well as the QP direct gaps for all the k-points
included in the calculation ({{kptgw}}).
To plot the difference QP - KS, use:

```{code-cell} ipython3
sigres.plot_qpgaps();
```

For the absolute QP gaps:

```{code-cell} ipython3
sigres.plot_qpgaps(plot_qpmks=False);
```

To plot the QP results as a function of the initial KS energy:

```{code-cell} ipython3
sigres.plot_qps_vs_e0(); #tight_layout=True);
```

and use `with_fields` to filter the quantity of interest:

```{code-cell} ipython3
sigres.plot_qps_vs_e0(with_fields=["vxcme", "sigxme", "sigcmee0"], sharey=True);
```

To plot the QP energies on top of the KS energies used in the SIGMA run:

```{code-cell} ipython3
sigres.plot_qpbands_ibz();
```

To plot the KS band structure with markers whose size is proportional to the QP correction
and whose direction gives the sign of the correction:

```{code-cell} ipython3
sigres.plot_ksbands_with_qpmarkers(fact=1000);
```

We can also plot the $<\Psi^{KS}_{mk}\,|\,\Psi^{QP}_{nk}>$ coefficients for given spin and k-point:

```{code-cell} ipython3
sigres.plot_eigvec_qp(spin=0, kpoint=0);
```

In this case, we have a diagonal matrix because the wavefunctions are not updated ($G_0W_0$).
The scenario is completely different if you start to perform self-consistent calculations with update
of the QP amplitudes.

+++

## Plotting the spectral function

This example shows how to plot the $G_0W_0$ spectral functions $A(\omega)$
at the $\Gamma$ point. See also lesson tgw2_4

```{code-cell} ipython3
with abilab.abiopen(abidata.ref_file("al_g0w0_sigmaw_SIGRES.nc")) as al_sigres:
    # Plot A(w) for the first spin, the gamma point, and bands in [0,1,2,3]
    al_sigres.plot_spectral_functions();
```

## Analyzing multiple SIGRES files with robots

To analyze the convergence of the QP results, we can use the SigresRobot.
Let's build our robot from a list of SIGRES files.

```{code-cell} ipython3
# List of SIGRES files computed with different values of nband.
filenames = [
    "si_g0w0ppm_nband10_SIGRES.nc",
    "si_g0w0ppm_nband20_SIGRES.nc",
    "si_g0w0ppm_nband30_SIGRES.nc",
]

filepaths = [abidata.ref_file(fname) for fname in filenames]

robot = abilab.SigresRobot.from_files(filepaths)
```

Then we plot the convergence of the QP direct gap as a function of the number of bands
in the self-energy for all the k-points available in the netcdf files:

```{code-cell} ipython3
robot.plot_qpgaps_convergence(sortby="sigma_nband", sharey=False);
```

If we are interested in the convergence of the real/imaginary part of the self-energy
and of the renormalization factor ...

```{code-cell} ipython3
robot.plot_qpdata_conv_skb(spin=0, kpoint=(0, 0, 0), band=3, sortby="sigma_nband");
```

We can also plot the QP data as a function of the KS energies on the same figure with:

```{code-cell} ipython3
robot.plot_qpfield_vs_e0("qpeme0", sortby="sigma_nband");
```

```{code-cell} ipython3
#robot.get_qpgaps_dataframe(spin=0, kpoint=(0, 0, 0))
```

<div class="alert alert-info" role="alert">
Robots can also be constructed from the command line with: abicomp.py sigres FILES
</div>
