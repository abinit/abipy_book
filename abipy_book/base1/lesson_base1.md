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

# Base1 lesson (H2 molecule)

<div class="jumbotron">
  <h1 class="display-3">First (basic) lesson with Abinit and AbiPy</h1>
  <p class="lead">The H<sub>2</sub> molecule</p>
  <hr class="my-4">
  <p>This lesson aims at showing how to get the following physical properties:</p>

   <ul>
    <li>the (pseudo) total energy</li>
    <li>the bond length</li>
    <li>the charge density</li>
    <li>the atomisation energy</li>
   </ul>
</div>

This tutorial is a complement to the standard [ABINIT tutorial on H$_2$](https://docs.abinit.org/tutorial/base1).
Here, powerful flow and visualisation procedures will be demonstrated.
Still, some basic understanding of the stand-alone working of ABINIT is a prerequisite.
Also, in order to fully benefit from this AbiPy tutorial, other more basic AbiPy tutorials should have been followed,
as suggested in the [introduction page](../intro).

There are three methodologies to compute the optimal distance between the two Hydrogen atoms.
One could:

   * compute the **total energy** for different values of the interatomic distance, make a fit through
     the different points, and determine the minimum of the fitting function;
   * compute the **forces** for different values of the interatomic distance, make a fit through
     the different values, and determine the zero of the fitting function;
   * use an automatic algorithm for minimizing the energy (or finding the zero of forces).

In this AbiPy notebook, we will be focusing on the first approach.
More specifically we will build an `Flow` to compute the energy and the forces in the $H_2$ molecule
for different values of the interatomic distance.
This exercise will allow us to learn how to generate multiple input files using AbiPy and
how to analyze multiple ground-state calculations with the AbiPy robots.

```{include} ../snippets/plotly_matplotlib_note.md
```

## Our first AbiPy function

```{code-cell}
import numpy as np

import warnings
warnings.filterwarnings("ignore")  # Ignore warnings

from abipy import abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

We need a function that generates an input file for a GS calculations for the $H_2$ molecule in a big box.
Ideally a function that receives the distance `x`, the cutoff energy `ecut` and the size of the big box
in input so that we can customize the output and generate multiple input objects easily.

Fortunately we already have such a function in the `lesson_base1.py` module.
Let's import it and look at the code:

```{code-cell}
from lesson_base1 import gs_input
abilab.print_source(gs_input)
```

If the function is called without arguments, the default values (specified in the prototype) are used.
Let's try:

```{code-cell}
gsinp = gs_input()
print("The value of ecut is:", gsinp["ecut"])
```

The `AbinitInput` is a dict-like object whose usage is documented in this [notebook](../abinit_input).
Inside jupyter, we can get the HTML representation of the input with:

```{code-cell}
gsinp
```

The input object can be converted into a string.
More importantly, an `AbinitInput` *has* an AbiPy structure (see [Structure notebook](../structure)),
a list of pseudopotential objects and provides several methods
to facilitate the specification of input variables.

```{code-cell}
print(gsinp.structure)
print("The big box volume is:", gsinp.structure.volume)
```

```{code-cell}
gsinp.structure.plot();
```

Let's print some info about our pseudopotentials:

```{code-cell}
print(gsinp.pseudos[0])
```

## Computation of the interatomic distance

+++

At this point, we can use `gs_input` to generate an [Abinit Flow](../flows)
to compute the total energy and the forces of H-H with different interatomic distances.
We have already prepared such a function in `build_flow`, let's have a look at the code:

```{code-cell}
from lesson_base1 import build_flow
abilab.print_source(build_flow)
```

Note that we are working at fixed {{ecut}} and {{acell}}, only the H-H distance is modified.
Let's call the function to build our flow:

```{code-cell}
flow = build_flow(options=None)
flow.get_graphviz()
```

Ok, so far so good.
With just three lines of codes and our `gs_input` function, we managed
to construct an AbiPy flow for the $H_2$ molecule.
Let's write some python code to check that we really obtained what we had in mind:

```{code-cell}
inputs = [task.input for task in flow.iflat_tasks()]

print("ecuts:\n", [inp["ecut"] for inp in inputs])

print("vols:\n", ["%.1f" % inp.structure.volume for inp in inputs])

def hh_dist(structure):
    return np.linalg.norm(structure.cart_coords[1] - structure.cart_coords[0])

from pprint import pprint
print("h-h [Ang]:\n", ["%.3f" % hh_dist(inp.structure) for inp in inputs])
```

At this point, we could run the flow in the notebook by just calling:

    flow.make_scheduler().start()

or, alternatively, execute the `lesson_base1.py` script to build
the directory with the flow and then use:

    abirun.py flow_h2 scheduler

inside the terminal.

Let's assume the flow has been already executed and let's focus on the analysis of the final results.

+++

## Analyzing the main output file

First of all, it is always a good idea to check whether the SCF cycle is converged.
Obviously one could open the main output file, find the SCF iterations and look for warnings but
there is a much faster (and better) way to do that with AbiPy:

```{code-cell}
abo = abilab.abiopen("flow_h2/w0/t0/run.abo")
print(abo)
```

To get the list of Warnings/Comments/Errors:

```{code-cell}
print(abo.events)
```

To plot the SCF cycle, use:

```{code-cell}
abo.plot();
```

Since this is not a structural relaxation, the initial and final structures must be equal:

```{code-cell}
abo.initial_structure == abo.final_structure
```

The basic dimensions and parameters of the run can be extracted from the output file with:

```{code-cell}
abo.get_dims_spginfo_dataset()
```

Within the shell, one can use:

    abiview.py abo flow_h2/w0/t0/run.abo

to plot the SCF cycle or

    abiopen.py flow_h2/w0/t0/run.abo

to open the file and start the ipython terminal

+++

Note that the object returned by `abiopen` has an associated file descriptor that should be closed
explicitly by the user with `abo.close()`.
This is the reason why many examples will use the python context manager:

```python
    with abiopen("path") as abifile:
        print(abifile)
```

that closes the file when we exit the `with` block.

+++

## Extracting results from the GSR files

The ground-state results are saved in the `GSR.nc` files whose API is extensively
discussed in the [GSR notebook](../gsr).

Let's have a look at the results produced by the first task:

```{code-cell}
with abilab.abiopen("flow_h2/w0/t0/outdata/out_GSR.nc") as gsr:
    print(gsr)
```

As we can see from the previous output, the `GSR` file contains information about
the crystalline structure, forces, stresses as well as the KS band structure.
In the jargon of object-oriented programming, one says that a `GsrFile` *has* a `Structure` object:

    gsr.structure

and *has* an `ElectronBands` object:

    gsr.ebands

This means that if you learn how to use the methods provided by `structure` and `ebands`, then you can
easily get these objects from the `GSR` file and use this API to post-process the results.
This is a general philosophy of AbiPy: every netcdf file object returned by `abiopen` contains
other objects (the structure is always available, while the presence of other objects depend of the particular file).
Remember this important point because we will use it a lot in the other lessons.

```{note}
To learn about the different methods provided by `structure` and `ebands`, you may want to
consult the API documentation

- {{Structure}}
- {{pymatgen_Structure}}
- {{ElectronBands}}
```

+++

Ok, now we know how to open and extract information from one `GSR` file.
In this tutorial, however, we need to analyze multiple `GSR` files!
If you are familiar with python, it should not be difficult to write a `for loop` that
iterates over a list of GSR files, extracts the total energy with the corresponding volume and creates two
lists that can be used to plot $E(d)$ with d the H-H distance.
This kind of operations are, however, very common and AbiPy provides a high-level interface (`robots`) to
operate on multiple files and post-process the data.

In the simplest case, one can use the `from_dir` method to tell the  `Robot` to
find all files of a particular type located within a directory tree,
and store all the data in memory:

```{code-cell}
robot = abilab.GsrRobot.from_dir("flow_h2")
robot
```

Once we have constructed the robot, we can start to invoke methods to extract/post-process the results.
For instance, we can construct a pandas Dataframe with:

```{code-cell}
table = robot.get_dataframe()
```

The table contains several columns:

```{code-cell}
table.keys()
```

Inside the notebook, we can visualize the table with:

```{code-cell}
table
```

Great! We managed to get a nice table with lot of useful results with just three lines of code and the robot!
There are however two problems:

   - The rows of the table are not ordered by increasing H-H distance (files are sorted alphabetically)

   - Our dataframe contains the energy of the different configurations, but we would like to plot the energy
     as a function of the H-H distance

Well, robots can do a lot of hard work but they are a little bit stupid so
we have to tell them what to do with the data.
More specifically we need a way to tell the robot that, for each `GSR` file, it should get the crystalline
structure, compute the distance between the first and the second atom and insert the result
in our table in a given column.
This kind of tasks are usually executed with `callbacks` i.e. functions that are passed in input
and **automatically executed** by the framework at runtime.

<img src="https://github.com/abinit/abipy_assets/blob/master/caution-machine_has_no_brain.jpg?raw=true" alt="">

Let's look at the documentation of `robot.get_dataframe`:

```{code-cell}
abilab.print_doc(robot.get_dataframe)
```

It seems complicated but the actual implementation of the callback is just three lines of code:

```{code-cell}
---
code_folding: []
run_control:
  marked: true
slideshow:
  slide_type: '-'
---
def hh_dist(gsr):
    """
    This callback receives a GSR file and computes the H-H distance.
    The robot will call this function to compute the H-H distance,
    and return a (key, value) tuple that will be inserted in the pandas DataFrame.
    """
    cart_coords = gsr.structure.cart_coords
    d = np.linalg.norm(cart_coords[1] - cart_coords[0])
    return "hh_dist", d

with abilab.GsrRobot.from_dir("flow_h2") as robot:
    table = robot.get_dataframe(funcs=hh_dist)
    table = table.sort_values(by="hh_dist")
```

As expected, now the table contains a new column with `hh_dist` in Angstrom:

```{code-cell}
"hh_dist" in table
```

Let's print the two columns with the H-H distance and the total energy:

```{code-cell}
table[["hh_dist", "energy"]]
```

Note that the energy in our `DataFrame` is given in eV to facilitate the integration
with `pymatgen` that uses eV for energies and Angstrom for lengths.
Let's add another column to our table with energies in Hartree:

```{code-cell}
table["energy_Ha"] = table["energy"] * abilab.units.eV_to_Ha
```

and use the `plot` method of pandas `DataFrames` to plot `energy_Ha` vs `hh_dist`

```{code-cell}
table.plot(x="hh_dist", y="energy_Ha", style="-o");
```

At this point, it should be clear that to plot the maximum of the forces as a function of the H-H distance
we just need:

```{code-cell}
table.plot(x="hh_dist", y="max_force", style="-o");
```

Want to plot the two quantities on the same figure?

```{code-cell}
table.plot(x="hh_dist", y=["energy_Ha", "max_force"], subplots=True);
```

```{note}

Your boss understands the data only if it is formatted inside a $\LaTeX$ tabular environment?

    table[["hh_dist", "energy"]].to_latex()

Need to send data to Windows users?

    table.to_excel("'output.xlsx")

Want to copy the dataframe to the system clipboard so that one can easily past the data
into an other applications e.g. Excel?

    table.to_clipboard()
```

## Analysis of the charge density

The `DEN.nc` file stores the density in real space on the FFT mesh.
A `DEN.nc` file *has* a `structure`, and an `ebands` object with the electronic eigenvalues/occupations
and a `Density` object with $n(r)$ (numpy array `.datar`) and $n(G)$ (`.datag`).

Let's open the file with `abiopen` and print it:

```{code-cell}
with abilab.abiopen("flow_h2/w0/t10/outdata/out_DEN.nc") as denfile:
    print(denfile)
    density = denfile.density
```

The simplest thing we can do now is to print $n(r)$ along a line passing through two points specified
either in terms of two vectors or two integers defining the site index in our `structure`.
Let's plot the density along the H-H bond by passing the index of the two atoms:

```{code-cell}
density.plot_line(0, 1);
```

Great! If we have a netcdf file and AbiPy, we don't need to use `cut3d` to extract the data from the file
and we can do simple plots with matplotlib.
Unfortunately, $n(r)$ is a 3D object and the notebook is not the most suitable tool to visualize this kind of dataset.
Fortunately there are several graphical applications to visualize 3D fields in crystalline environments
and AbiPy provides tools to export the data from netcdf to the text format supported by the external graphical tool.

For example, one can use:

```{code-cell}
#density.visualize("vesta")
```

to visualize density isosurfaces of our system:

![](https://github.com/abinit/abipy_assets/blob/master/h2_density.png?raw=true)

+++

## Conclusions

To summarize, we learned how to define python functions that can be used to generate many input files easily.
We briefly discussed how to use these inputs to build a basic AbiPy flow without dependencies.
More importantly, we showed that AbiPy provides several tools that can be used to inspect and analyze
the results without having to pass necessarily through the creation and execution of the `Flow`.
Last but not least, we discussed how to use `robots` to collect results from the output files and store
them in pandas DataFrames.

AbiPy users are **strongly recommended** to familiarize themself with this kind of interface before
moving to more advanced features such as the flow execution that requires a good understanding of the python language.
As a matter of fact, we decided to write AbiPy in python not for efficiency reasons (actually python
is usually slower that Fortran/C) but because there are tons of libraries for scientific applications
(numpy, scipy, pandas, matplotlib, jupyter, etc).
If you learn to use these great libraries for your work you can really boost your productivity and save a lot of time.

+++

A logical next lesson would be the the tutorial about the [ground-state properties of silicon](../base3/lesson_base3)
