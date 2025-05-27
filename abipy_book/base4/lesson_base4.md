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

# Base4 lesson (Aluminium)

<div class="jumbotron">
  <h1 class="display-3">Fourth (basic) lesson with Abinit and AbiPy</h1>
  <p class="lead">The H<sub>2</sub> molecule</p>
  <hr class="my-4">
  <p>This lesson aims at showing how to get the following physical properties, for a metal, and for a surface:

the total energy
the lattice parameter
the relaxation of surface atoms
the surface energy You will learn about the smearing of the Brillouin zone integration, and also a bit about preconditioning the SCF cycle.
This lesson should take about 1 hour and 30 minutes.
</p>
  <p class="lead">
    <a class="btn btn-primary btn-lg" href="#" role="button">Learn more</a>
  </p>
</div>This tutorial is a complement to the standard [ABINIT tutorial on aluminum](https://docs.abinit.org/tutorial/base4). 
Here, powerful flow and visualisation procedures
will be demonstrated. Still, some basic understanding of the stand-alone working of ABINIT is a prerequisite.
Also, in order to fully benefit from this Abipy tutorial, other more basic Abipy tutorials should have been followed,
as suggested in the [abitutorials index page](https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb).

```{include} ../snippets/plotly_matplotlib_note.md
```

```{include} ../snippets/manager_note.md
```

```{code-cell}
import warnings
warnings.filterwarnings("ignore")  # Ignore warnings

import numpy as np

from abipy import abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook

# This line tells the notebook to show plots inside of the notebook
%matplotlib notebook
```

## The convergence study with respect to k points and broadening

Note that there is usually a STRONG cross-convergence effect between the number
of k-points and the value of the broadening, tsmear.
The right procedure is: for each value of tsmear, reach convergence with respect
to the number of k points, then compare the k-point converged values for different values of tsmear.

In what follows, we will restrict ourselves to the grids with nkpt=2, 10 and 28.

+++

As usual, we start by writing an helper function to generate the input for the structural relaxation of Al.
The function accepts `tsmear` and another parameter, `nksmall`, that will be used to define the BZ sampling:

```{code-cell}
from lesson_base4 import relax_input
abilab.print_source(relax_input)
```

Now we use `relax_input` to generate multiple inputs with different values of `tsmear` and `nksmall`
and we pass the input objects to the Flow constructor.
To keep things as simple as possible, we use independent tasks:

```{code-cell}
from lesson_base4 import build_relax_tsmear_nkpts_convflow
abilab.print_source(build_relax_tsmear_nkpts_convflow)
```

```{code-cell}
flow = build_relax_tsmear_nkpts_convflow(options=None)
flow.get_graphviz()
```

```{code-cell}
for task in flow.iflat_tasks():
    print(task.pos_str, "ngkpt", task.input["ngkpt"],
          "tsmear", task.input["tsmear"], "occopt", task.input["occopt"])
```

```{code-cell}
#abo = abilab.abiopen("flow_al_relax/w0/t0/run.abo")
#print(abo)
```

```{code-cell}
#abo.plot();
```

```{code-cell}
#hist = abilab.abiopen("flow_al_relax/w0/t0/outdata/out_HIST.nc")
#print(hist)
```

```{code-cell}
#hist.plot();
```

```{code-cell}
robot = abilab.GsrRobot.from_dir("flow_al_relax_tsmear_nkpt")
robot
```

```{code-cell}
robot.remap_labels(lambda gsr: "nkpt: %d, tsmear %.2f" % (gsr.nkpt, gsr.tsmear))
```

```{code-cell}
data = robot.get_dataframe()
print(data.keys())
```

First of all, let's sort our rows first by `nkpt` and then by `tsmear` inside each `nkpt` group so that
we have print the table in a nice format:

```{code-cell}
data = data[["nkpt", "tsmear", "a", "alpha", "energy", "pressure", "volume", "max_force"]]
```

```{code-cell}
data = data.sort_values(by=["nkpt", "tsmear"])
data
```

```{code-cell}
import seaborn as sns
sns.pairplot(data, x_vars="nkpt", y_vars=["energy", "a", "volume"], hue="tsmear");
```

```{code-cell}
y_vars = ["energy", "structure.lattice.a", "structure.volume"]

robot.plot_convergence_items(y_vars, sortby="nkpt", hue="tsmear");
```

```{code-cell}
eframe = robot.get_energyterms_dataframe();
eframe
```

```{code-cell}
with abilab.abiopen("flow_al_relax_tsmear_nkpt/w0/t11/outdata/out_GSR.nc") as gsr:
    ebands_w0t11 = gsr.ebands

ebands_w0t11.plot();
```

```{code-cell}
ebands_w0t11.boxplot();
```

Now you might ask yourself: "The total energy with nkpt == 2 is clearly not converged wrt tsmear. 
What are the effects of the smearing on the KS eigenvalues for nkpt == 2?"
The `GsrRobot` can construct an `ElectronBandsPlotter` that allows us to compare multiple band structures
so it's just a matter of telling the robot that we want a plotter object in which only the
`GSR` files with only two k-points in the IBZ, then we can use the plotter to visualize the results:

```{code-cell}
plotter = robot.get_ebands_plotter()
plotter.gridplot();
```

```{code-cell}
robot.gridplot_with_hue(hue="nkpt");
```

The plotter provides different plot methods to visualize the same data.
Perhaps you prefer this:

```{code-cell}
plotter.combiboxplot();
```

```{code-cell}
#robot.boxplot_ebands();
```

```{code-cell}
robot.combiboxplot_ebands();
```

```{code-cell}
edos_plotter = robot.get_edos_plotter(width=0.1, step=0.2)
```

```{code-cell}
edos_plotter.combiplot();
```

```{code-cell}
edos_plotter.gridplot();
```

```{code-cell}
#e3d = ebands_w0t11.get_ebands3d()
#e3d.plot_isosurfaces()
```

A logical next lesson would be the the tutorial on
[phonon calculations with DFPT](https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/dfpt/lesson_dfpt.ipynb)

+++
