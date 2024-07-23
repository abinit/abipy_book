---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python
  language: python
  name: python3
---

# Bethe-Salpeter calculations with AbiPy

This lesson discusses how to calculate the macroscopic dielectric function, $\epsilon_\infty(\omega)$,
including excitonic effects within the Bethe-Salpeter equation (BSE).
Crystalline silicon is used as test case.

For a more detailed description of the Abinit implementation, see the official Abinit
[BSE tutorial](https://docs.abinit.org/tutorial/bse/).
A brief description of the formalism can be found in the [BSE_notes](https://docs.abinit.org/theory/bse/).

## Typical BSE flowchart

The flowchart of a typical Bethe-Salpeter run is schematically depicted in the diagram below:

+++

<img src="bse_flowchart.svg" alt="bse_flowchart">

+++

The WFK file (KSS file in old versions of Abinit) contains the Kohn-Sham (KS) wavefunctions
and energies and is represented with an ellipsis.
The path on the left indicated with blue arrows represents the RPA calculation ({{optdriver}} = 3)
that produces the SCR file (see also the first lesson of the $GW$ tutorial).
Once the WFK (KSS) and the SCR file are available, we can finally contruct the BSE Hamiltonian
and solve the Bethe-Salpeter problem (the green rectangle at the bottom of the flowchart).
The construction of the Bethe-Salpeter Hamiltonian represents a significant portion
of the overall CPU time due to the large number of transitions (bands and in particular $k$-points)
needed for an accurate description of the frequency-dependence of the polarizability.

For BSE computations, it is common practice to simulate the self-energy corrections
by employing the scissors operator whose value can be obtained either from experiments
or from *ab-initio* calculations.
The scissors operator allows one to avoid a costly $GW$ calculation that should performed
for all the $k$-points and bands included in the transition space
(the optional path on the right indicated with yellow arrows that corresponds to ({{optdriver}} = 4).

For this reason, in this lesson, we will employ two commonly used approximations that will
reduce considerably the computational cost of the BSE flowchart while giving reasonably accurate results:

   * The *ab-initio* $W$ is replaced by a model dielectric function
     that is constructed from the GS density $n(r)$ and the additional variable {{mdf_epsinf}}
     that gives the value of the static limit $\epsilon_\infty(\omega=0)$.
     This approximation allows us to bypass the blue boxes in the diagram above ({{optdriver}} = 3)

   * The modifications introduced by the $GW$ self-energy on the initial KS band structure
     are approximated with a scissor operator {{mbpt_sciss}}.
     This approximation allows us to bypass the yellow boxes in the diagram above (`optdriver=4`).

Under these assumptions, the BSE flowchart reduces to a simple GS-SCF run to get $n(r)$ plus
a NSCF calculation of the band structure on a dense $k$-mesh and, finally, the solution
of the BSE problem (the green box).

+++

## AbiPy flow for BSE with the model dielectric function

After the standard imports:

```{code-cell}
import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

import abipy.abilab as abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook
import abipy.flowtk as flowtk

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

We import from `lesson_bse` the function that builds the 3 input objects we are going
to use to build the `Flow`:

```{code-cell}
from lesson_bse import make_scf_nscf_bse_inputs
abilab.print_source(make_scf_nscf_bse_inputs)
```

```{code-cell}
scf_inp, nscf_inp, bse_inp = make_scf_nscf_bse_inputs(ngkpt=(4, 4, 4), ecut=6, ecuteps=3)
```

The function `make_scf_nscf_bse_inputs` returns three `AbinitInput` objects:

   * The first input (`scf_inp`) solves the KS equations on a Monkhorst-Pack mesh
     to obtain the groud-state density $n(r)$.

   * The second input (`nscf_inp`) uses the density produced by `scf_inp` to compute
     the KS band structure on a **randomly-shifted** $k$-mesh in order to accelerate
     the convergence of the optical properties with respect to the $k$-sampling.

   * Finally, the third input (`bse_inp`) uses the `WFK` file produced in the previous step
     to solve an approximated BSE equation in which the ab-initio
     screened interation $W$ is approximated
     by a model dielectric function that depends only on $n(r)$ and the input variable {{mdf_epsinf}}
     that gives the value of $\epsilon_\infty(0)$

+++

The variables governing the BSE run are those in the `vargw` section of `bse_inp`:

```{code-cell}
bse_inp
```

Once we have our three input objects, we can create a flow to automate the calculation.
Note that AbiPy already provides the `BseMdfWork` class that is explicitly designed for this kind of calculation:

```{code-cell}
from lesson_bse import build_bse_flow
abilab.print_source(build_bse_flow)
```

Let's build the flow:

```{code-cell}
flow = build_bse_flow(options=None)
```

The graphical representation of the flow reveals that the `BseTask` depends on the `NscfTask` that
in turns depends on the initial `ScfTask`.

```{code-cell}
flow.get_graphviz()
```

```{code-cell}
#flow.plot_networkx(with_edge_labels=True);
```

If you are working with python, you can build the directories of the Flow with:

    flow.build_and_pickle_dump()

+++

Now you can execute the *lesson_bse.py* script to generate the flow  and then use:

    abirun.py flow_bse scheduler

```{include} ../snippets/abicheck_warning.md
```

Alternatively, one can use the files in the github repository and use AbiPy
to analyze the data.

+++

## Analyzing the results

Now we can finally analyze the results. In this case, we are mainly interested in the
frequency-dependent macroscopic dielectric function, $\epsilon_\infty(\omega)$, produced
by the `BseTask`

```{code-cell}
# The BseTask is the last task in the first work
# i.e. flow[0][2] or, much better, flow[0][-1]
bse_task = flow[0][-1]
bse_task
```

The `BseTask` has produced a netcdf file (`MDF.nc`) containing the most important results of the run.
Let's open the file with:

```{code-cell}
---
run_control:
  marked: true
---
mdf_file = abilab.abiopen("flow_bse/w0/t2/outdata/out_MDF.nc")
print(mdf_file)
```

and use `matplotlib` to plot the imaginary part of $\epsilon_\infty(\omega)$

```{code-cell}
mdf_file.plot_mdfs();
```

Meaning of the three curves:

   * EXC is $\epsilon_\infty(\omega)$ computed from the BSE with excitonic effects included
   * KS-RPA is the analogous quantity computed within the RPA and the KS band structure
   * GW-RPA corresponds to the RPA expression but computed with modified band energies obtained
     by "opening" the KS eigenvalues with a constant scissor operator that tries to mimic
     the $GW$ corrections (`soenergy` variable)

It is worth stressing that:

1) The RPA-KS spectrum underestimates the experimental optical threshold
   due to the well-know band-gap problem of DFT.
   Most importantly, the amplitude of the first peak is underestimated.

2) The RPA-GW results with QP corrections simulated with `soenergy` does not show
   any significant improvement over RPA-KS: the RPA-GW spectrum is just shifted towards
   higher frequencies due to opening of the gap, but the shape of the two spectra
   is very similar, in particular the amplitude of the first peak is still underestimated.

3) On the contrary, the inclusion of the BSE kernel leads to important changes both
   in the optical threshold as well as in the amplitude of the first peak.
   This simple analysis tells us that the first peak in the absorption spectrum of silicon
   has a strong excitonic character that is not correctly described within the RPA.
   Our first BS spectrum is not converged at all and it barely resembles
   the experimental result, nevertheless this unconverged calculation is already able
   to capture the most important physics.

The difference among the three approaches is schematically depicted in the figure below:

+++

<img src="schematic_dft_gw_bse.svg" alt="schematic dft gw bse">

+++

To plot the real part of $\epsilon_\infty(\omega)$

```{code-cell}
mdf_file.plot_mdfs(cplx_mode="re");
```

It should be stressed that the screened interaction $W$ is the fundamental ingredient that leads
to the attractive interaction between electrons and holes (excitonic effects).
In a metallic system, the dielectric function is large, $W$ is small and excitonic effects are strongly damped.

To understand this point, we can do a test calculation with a very large value of {{mdf_epsinf}}
so that our BSE Hamiltonian will be constructed with a metallic $W$:

```{code-cell}
from lesson_bse import build_bse_metallicW_flow
abilab.print_source(build_bse_metallicW_flow)
```

```{code-cell}
metalW_flow = build_bse_metallicW_flow(options=None)
```

Let's assume we have already executed the flow and let's have a look at the results:

```{code-cell}
with abilab.abiopen("flow_bse_metallicW/w0/t2/outdata/out_MDF.nc") as mdf_file:
    mdf_file.plot_mdfs();
```

As you can see, the EXC curve computed with a metallic $W$
is similar to the results obtained in GW-RPA.
In particular the first peak is now shifted towards higher frequencies and its amplitude
is decreased when compared to the previous results.

This behaviour can be easily understood if we consider that
the BSE formalism reduces to the RPA if $W$ tends to 0.
The EXC curve is still shifted towards higher frequencies when compared with KS-RPA
but this effect is mainly due to the scissor operator that opens the KS gap.
A similar calculation done with `soenergy=0` would give an EXC curve similar to KS-RPA.
This test is left as an optional exercise.

+++

## Convergence study with respect to the $k$-point sampling

The most important parameter that should be checked for convergence is the number of $k$-points.
This convergence study represents the most tedious and difficult part since it requires
the generation of new WFK files for each k-mesh
(the list of $k$-points for the wavefunctions and the set of $q$-points in the screening
must be consistent with each other).

In the previous section, we have shown how to build a flow for BSE calculation with a fixed
$k$-points sampling.
We can thus reuse the same logic to construct a `Flow` made of multiple `BseMdfWorks`,
each `Work` will have a different $k$-point sampling.

Let's create, for example, a `Flow` that solves that BSE equation on
a `4x4x4`, `6x6x6` and a `8x8x8` $k$-mesh:

```{code-cell}
from lesson_bse import build_bse_kconv_flow
abilab.print_source(build_bse_kconv_flow)
```

```{code-cell}
flow_kconv = build_bse_kconv_flow(options=None)
```

```{code-cell}
flow_kconv.get_graphviz()
```

```{code-cell}
#flow_kconv.plot_networkx();
```

Change the *lesson_bse.py* script so that `build_bse_kconv_flow` is called in main instead of `build_bse_flow`.
Run the script and submit the calculation with abirun.py FLOWDIR scheduler as usual.

+++

Our `Flow` has three `BseTasks` and therefore three different `MDF.nc` files containing $\epsilon_\infty(\omega)$.
The MDF files are available in the github repository.

In order to plot the three $\epsilon_\infty(\omega)$ on the same graph, we have use the `MdfRobot`
that will gather the results for us:

```{code-cell}
robot = abilab.MdfRobot.from_dir("flow_bse_kconv")
robot
```

We can now finally compare the imaginary and the real part of $\epsilon_\infty(\omega)$ with `matplotlib`:

```{code-cell}
plotter = robot.get_multimdf_plotter()
```

```{code-cell}
plotter.plot();
```

## Exercises

Use `make_scf_nscf_bse_inputs` and `BseMdfWork` to perform the following convergence studies:

   * Convergence with respect to the number of planewaves in the screening {{ecuteps}}
   * Convergence with respect to the number of $k$-points

See also the discussion reported in the official [BSE tutorial](https://docs.abinit.org/tutorial/bse/)
