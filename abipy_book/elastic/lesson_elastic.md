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

# Elastic and piezoelectric properties with Abinit and AbiPy

+++

This tutorial shows how to calculate the following physical properties related to strain:

  * the *clamped-ion* elastic tensor
  * the *clamped-ion* piezoelectric tensor (insulators only)
  * the *internal* strain tensor
  * the atomic relaxation corrections to the elastic and piezoelectric tensor (*relaxed-ion tensors*)

The discussion is based on the [tutorial on elastic properties](https://docs.abinit.org/tutorial/elastic/)
available on the Abinit web site.
More specifically, we will discuss how to

   * build a flow to compute all the ingredients needed for elastic and piezoelectric tensors
   * use anaddb and AbiPy to obtain the tensors and compute several important physical properties
     starting from the final DDB file produced by the flow
   * perform a convergence study with respect to the k-mesh and use the `DdbRobot` to analyze the convergence.

You might find additional material related to the present section
in the following references (already mentioned in the official tutorial):

* [Systematic treatment of displacements, strains, and electric fields in density-functional perturbation theory](https://dx.doi.org/10.1103/physrevb.72.035105)
* [Metric tensor formulation of strain in density-functional perturbation theory](https://doi.org/10.1103/physrevb.71.035117)

The first paper provides a detailed discussion of the theory underlying the incorporation
of atom-relaxation corrections.
We strongly recommend to read this article as this notebook will mainly focus on the usage
of the python API assuming you are already familiar with the theoretical aspects.
The second paper discusses in more details the DFPT treatment of strain perturbations in Abinit.

If you are already familiar with python and AbiPy-Abinit are already installed and configured,
you may want to use directly the command line interface.
There is a README.md file in the directory of this lesson explaining how to analyze the data from the shell
using ipython and matplotlib.

```{note}
The code in this notebook requires abinit >= 8.9 and abipy >= 0.6
```

## DFPT calculation of elastic and piezolectric tensors

Before starting, we need to import the python modules and the functions needed in the notebook:

```{code-cell}
import numpy as np
import warnings
warnings.filterwarnings("ignore") # to get rid of deprecation warnings

from abipy import abilab
abilab.enable_notebook() # This line tells AbiPy we are running inside a notebook
import abipy.flowtk as flowtk

# This line configures matplotlib to show figures embedded in the notebook.
# Replace `inline` with `notebook` in classic notebook
%matplotlib inline

# Option available in jupyterlab. See https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget
```

and an useful function from the `lesson_elastic` module required to generate our DFPT flows:

```{code-cell}
from lesson_elastic import make_scf_input
abilab.print_source(make_scf_input)
```

The function makes some assumptions for important parameters such as
the crystalline structure and the pseudos.
This is done on purpose to keep the code as simple as possible.
It should not be so difficult to generalize the implementation to take into account other cases.
Note how the function accepts an optional argument {{ngkpt}} defining the **k**-mesh so that we can easily
change the sampling *e.g.* for convergence studies.

Let's start to play with our new function:

```{code-cell}
scf_input = make_scf_input()
scf_input
```

```{code-cell}
print(scf_input.structure)
```

```{code-cell}
scf_input.structure.plot();
```

We are using the same norm-conserving pseudopotentials of the official tutorial but
this does not mean you should use them for production calculations (there must be a reason
why the directory is called **Psps_for_tests**).

There are 16 valence electrons per unit cell hence {{nband}} has been set to 8
(yes, we are dealing with a non-magnetic semiconductor):

```{code-cell}
for pseudo in scf_input.pseudos:
    print(pseudo, "\n")
```

The `scf_input` represents the **building block** for our DFPT calculation.
As usual, AbiPy provides a `Work` to compute elastic and piezoelectric properties
starting from an input representing a ground-state calculation
thus it is just a matter of calling `make_scf_input` and pass the result to
`ElasticWork.from_scf_input` to construct our flow.

Let's have a look at the actual implementation:

```{code-cell}
from lesson_elastic import build_flow
abilab.print_source(build_flow)
```

Now we can call the function to build our flow:

```{code-cell}
flow = build_flow()
flow.show_info()
```

and use the `get_graphviz` method to visualize the connection among the `Tasks`:

```{code-cell}
flow.get_graphviz()
```

```{code-cell}
# matplotlib version based on networkx
#flow.plot_networkx(with_edge_labels=True);
```

In a nutshell:

   * we compute the `WFK` file in the `ScfTask` (red circle)
   * the ground-state wavefunctions are used by the three `DdkTasks` to compute
     $\dfrac{\partial u}{\partial{\bf k}}$ for the three different directions.
   * The `ElasticTasks` needs the `WFK` file to compute the six strain perturbations
     (3 for uniaxial and 3 for shear strain) while the `DDK` files are required
     to compute the mixed 2nd-order derivatives with respect to strain and electric field
     needed for the piezoelectric tensor.

Note that, **contrarily to the approach used in the standard tutorial**, the AbiPy `Work` does not use datasets.
The perturbations of interest (strain, atomic-perturbation, ddk, electric field)
are obtained with different `Tasks` that can be executed in parallel.

To understand better this point, we can print a table with the most important Abinit variables defining
the DFPT calculation.
If the variable is not present in the input, the entry is set to `None`.

```{code-cell}
flow.get_vars_dataframe("qpt", "rfphon", "rfatpol", "rfdir", "rfelfd", "rfstrs", "kptopt")
```

If the meaning of these variables is not clear, you can consult the [Abinit documentation](https://docs.abinit.org)
or access the documentation directly from python with *e.g.*:

```{code-cell}
abilab.docvar("rfstrs")
```

```{note}
For your convenience the links to the doc of the different variables are listed below:

- {{qpt}}
- {{rfphon}}
- {{rfatpol}}
- {{rfdir}}
- {{rfelfd}}
- {{rfstrs}}
- {{kptopt}}
```

Now we can generate the `flow_elastic` directory with the input files by executing the `lesson_elastic.py` script.
Then use the {{abirun}} script to launch the entire calculation with:

    abirun.py flow_elastic scheduler

You will see that all `PhononTasks` and `ElasticTasks` are executed in parallel on your machine
once the three `DdkTasks` are completed.

```{include} ../snippets/abicheck_warning.md
```

If you prefer to skip this part, you may want to jump to next section about the post-processing of the results.
Note that the output files are already available in the repository so it is also possible to try
the AbiPy post-processing tools without having to run the flow.

+++

## Post-processing the results

Our flow is completed and we have the final DDB file
in the `outdata` directory of the `Work` (AbiPy has automatically merged all the partial DDB files
at the end of the calculation by invoking `anaddb` for you).
Let's open this DDB file with:

```{code-cell}
ddb = abilab.abiopen("flow_elastic/w0/outdata/out_DDB")
print(ddb)
```

The `DdbFile` object provides an easy-to-use interface that invokes `anaddb` to post-process
the data stored in the DDB file.
All the methods that invoke `anaddb` use the `ana` prefix so it is not strange to
see that we can obtain the elastic and piezoelectric tensors by just calling:

```{code-cell}
edata = ddb.anaget_elastic(verbose=1)
```

Note that we are calling the method without arguments. This means that AbiPy will try to detect
**automatically** how to set the anaddb input variables.
The docstring of the method explains the logic used to set the variables.

```{code-cell}
abilab.print_doc(ddb.anaget_elastic)
```

Let's print the object to get a summary of the most important results:

```{code-cell}
print(edata)
```

Since the DDB file contains `internal strain terms` and piezoelectric terms,
AbiPy set `elaflag` to 3, `instrflag` to 1 and `piezoflag` to 3 so that anaddb
will compute both `relaxed` and `clamped-ion` tensors.

```{code-cell}
abilab.docvar("elaflag", executable="anaddb")
```

The tensors are available as attributes of the `edata` object:

```{code-cell}
edata.elastic_relaxed
```

These are pymatgen [tensors](https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/elasticity/tensors.py),
more specifically [ElasticTensor objects](https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/elasticity/elastic.py)
so we have access to several useful methods.
To get the Voigt bulk modulus, use:

```{code-cell}
edata.elastic_relaxed.k_voigt
```

while the compliance tensor is easily obtained with:

```{code-cell}
edata.elastic_relaxed.compliance_tensor
```

One can also use the [elate](http://progs.coudert.name/elate) online tool to analyse the elastic tensor.

+++

To build a pandas DataFrame with properties derived from the elastic tensor
and the associated structure, use:

```{code-cell}
edata.get_elastic_properties_dataframe(properties_as_index=True)
```

For the meaning of the different quantities please consult the
[pymatgen module](https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/elasticity/elastic.py).

+++

To construct a dataframe with the Voigt indices and the tensor elements use:

```{code-cell}
edata.get_elastic_voigt_dataframe(tol=1e-5)
```

Note that the Voigt indices are given following the Python (C) notation
in which we start to count from zero.

This might be a bit confusing, especially when comparing with results reported in the literature.
The reason why we opted with the 0-based notation is that it facilitates the integration between
the DataFrame and other python methods.
A similar approach is used in AbiPy when e.g. one has to specify the band or the phonon mode index.

+++

At this point, it should be clear how to analyze the *relaxed-ion* piezoelectric tensor:

```{code-cell}
edata.piezo_relaxed
```

```{code-cell}
edata.get_piezo_voigt_dataframe(tol=1e-6)
```

## Convergence study wrt the number of k-points

In this part of the tutorial, we discuss how to compute elastic and piezoelectric
tensors with different k-point meshes and how to use the `DdbRobot` to analyze the results.

In principle, one should relax the structure with different k-point samplings and then
use the relaxed structure to compute elastic and piezoelectric properties.
This is easy with AbiPy but, for the time being, we ignore this point and use
the same structure so that we can learn how to use Python to analyze multiple calculations.
At the end of this tutorial you will find an example of flow in which the structural relaxation
is performed with different k-meshes.

Since we do not have to change the structure, performing a convergence study with respect to k-points
is just a matter of creating multiple `Works` inside a loop over k-meshes (our `make_scf_input`
is already accepting `ngkpt` in input):

```{code-cell}
from lesson_elastic import build_ngkpt_convflow
abilab.print_source(build_ngkpt_convflow)
```

```{code-cell}
ngkpt_flow = build_ngkpt_convflow(options=None, ngkpt_list=([2, 2, 2], [4, 4, 4], [8, 8, 8]))
```

```{code-cell}
ngkpt_flow.get_graphviz()
```

```{code-cell}
#ngkpt_flow.get_vars_dataframe("ngkpt")
```

To generate the flow with the `lesson_elastic.py` script, open the file,
comment the call to `build_flow` and uncomment `build_ngkpt_convflow`.
Then run the script and launch the calculation with:

    abirun.py flow_elastic_ngkpt_conv scheduler

as usual.

There are several output files located inside the `outdata` directories:

```{code-cell}
!find flow_elastic_ngkpt_conv/ -name "*_DDB"
```

Remember that our goal is to analyze the convergence of the elastic and piezoelectric properties
as function of `nkpt`.
So we are mainly interested in the final DDB files located in the `outdata` directories
of the works (`w0/outdata`, `w1/outdata`, `w2/outdata`).
These are indeed the DDB files with all the information needed in anaddb to compute the tensors.

The code below tells our robot that we would like to analyze all the DDB files
located in the output directories of the works:

```{code-cell}
robot = abilab.DdbRobot.from_dir_glob("./flow_elastic_ngkpt_conv/w*/outdata/")
robot
```

The DDB file are available in `robot.abifile`.
Each `DdbFile` object has a header (dictionary) containing metadata extracted from the DDB file.
To get the keywords in the header, use:

```{code-cell}
robot.abifiles[0].header.keys()
```

We will be using these meta variables to construct our pandas Dataframe so that we can analyze
the convergence of our physical quantities with e.g. `nkpt`.

Let's call `anacompare_elastic` to construct a DataFrame (`data`) with the elastic properties obtained
with the three DDB files and add the value of `ddb.header["nkpt"]`.
`elastdata_list` is a list of `ElastData` object we can use afterwards to access the individual tensors:

```{code-cell}
data, elastdata_list = robot.anacompare_elastic(ddb_header_keys="nkpt")
```

```{code-cell}
data.keys()
```

```{code-cell}
data
```

```{code-cell}
data[["k_voigt", "nkpt", "tensor_name"]]
```

To plot the convergence of selected properties versus the number of k-points, use:

```{code-cell}
robot.plot_xy_with_hue(data, x="nkpt", y=["k_voigt", "g_voigt", "y_mod"], hue="tensor_name");
```

For more examples on the use of DDB and robots, see the
[DDB notebook](https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/ddb.ipynb)

Now let's try something a bit more complicated.
Assume we want to plot the convergence of the individual elements of the tensor as as function of `nkpt`.
In this case, we have to work a bit more to create a Dataframe with the elements in Voigt notation,
add the value of `nkpt` associated to this tensor and finally concatenate the results in a single DataFrame:

```{code-cell}
df_list = []
for edata, ddb in zip(elastdata_list, robot.abifiles):

    # Get dataframe with tensor elements in Voigt notation
    df = edata.get_elastic_voigt_dataframe()

    # Add metadata and store dataframe in df_list
    df["nkpt"] = ddb.header["nkpt"]
    df_list.append(df)

# Concatenate dataframes
import pandas as pd
data = pd.concat(df_list)
data.head()
```

To select only the (0, 0) elements, use the syntax:

```{code-cell}
c00 = data[data["voigt_cinds"] == (0, 0)]
c00
```

Now we can finally plot the (0, 0) tensor elements as function of `nkpt` with:

```{code-cell}
c00.plot(x="nkpt", y=[k for k in c00 if k.startswith("elastic_")], subplots=True, style="-o");
```
