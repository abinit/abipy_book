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

# Phonons and Born effective charges

This lesson discusses how to compute phonon band structures, DOS 
and Born effective charges with Abinit and AbiPy.
The discussion closely follows the [second lesson](https://docs.abinit.org/tutorial/rf2/index.html)
on DFPT available on the Abinit web site. 
More specifically, we will discuss how to

   * Perform a convergence study for the phonon frequencies at $\Gamma$ as function of `ecut`
   * Compute the full phonon band structure of `AlAs` with the inclusion of LO-TO splitting
   * Obtain thermodynamic properties within the harmonic approximation

We assume that you have read the references mentioned in the [first Abinit lesson](https://docs.abinit.org/tutorial/rf1/index.html)
on DFPT.
You might find additional material, related to the present section, in the following references: 
   
* [Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.10355)
* [Phonons and related crystal properties from density-functional perturbation theory](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.73.515)

If you are already familiar with python and AbiPy-Abinit are already installed and configured,
you may want to use directly the command line interface.
See the README.md file in the directory of this lesson explaining how to analyze the data from the shell
using ipython and matplotlib.

## Phonon frequencies at $\Gamma$ as function of ecut

Before starting, we need to import the python modules and the functions we will need in the notebook:

```{code-cell} ipython3
# Use this at the beginning of your script so that your code will be compatible with python3
from __future__ import print_function, division, unicode_literals

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

and an useful function from the `lesson_dfpt` module that will be used to generate our DFPT flows:

```{code-cell} ipython3
from lesson_dfpt import make_scf_input
abilab.print_source(make_scf_input)
```

The function makes some assumptions for important parameters such as 
the crystalline structure and the pseudos. 
This is done on purpose to keep the code as simple as possible.
It should not be so difficult to generalize the implementation to take into account other cases.
Let's start to play with our new function:

```{code-cell} ipython3
scf_input = make_scf_input()
scf_input
```

```{code-cell} ipython3
print(scf_input.structure)
```

```{code-cell} ipython3
scf_input.structure.plot();
```

We are using the same pseudopotentials than the official tutorial.
Note that the xc functional is LDA in both pseudos but with a different
parametrization. This is the reason why we are using `ixc` in the input file.

```{code-cell} ipython3
for pseudo in scf_input.pseudos:
    print(pseudo, "\n")
```

As you can see, we essentialy have a standard input to perform a GS calculation. This object will represent
the **building block** for our DFPT calculation with AbiPy.

It this is not your first time you use the DFPT part of Abinit, you already know that phonon calculations
require an initial GS run to produce the `WFK` file 
followed by a DFPT run that reads the `WFK` file and solves the Sternheimer equations for $N_{\text{irred}}(q)$ 
atomic perturbations where $N_{\text{irred}}$ is the number of independent atomic displacements (assuming $q$ belongs to the k-mesh).

If you try to do a convergence study wrt `ecut` **without multi-datasets**, you will likely start from an initial GS input file with a given value of `ecut`, use it as a template to generate the DFPT input files, create symbolic 
links to the `WFK` file produced in the first GS step and then instruct Abinit to read this file with `irdwfk`.
Once you have a set of input files that work for a particular `ecut`, one can simply replicate the set of 
directories and files and use a script to change the value of `ecut` in the input files.
Then, of course, one has to run the calculations manually, collect the results and produce nice plots to understand
what is happening.

This approach is obviously boring and error-prone if you are a human being but it is easy to implement in an algorithm 
and machines do not complain if they have a lot of repetive work to do!
There are also several **technical advantages** in using this **task-based approach vs multi-datasets** but we discuss this point in more details afterwards. 

If the machine could speak, it will tell you: give me an object that represents an input for GS calculations,
give me the list of q-points you want to compute as well as the parameters that must be changed in the initial input 
and I will generate a `Flow` for DFPT calculations.
This logic appears so frequenty that we decided to encapsulate it in the `flowtk.phonon_conv_flow` factory function: 

```{code-cell} ipython3
from lesson_dfpt import build_flow_alas_ecut_conv
abilab.print_source(build_flow_alas_ecut_conv)
```

Let's call the function to build our flow:

```{code-cell} ipython3
flow = build_flow_alas_ecut_conv(options=None)

flow.show_info()
```

and call the `get_graphviz` method to visualize the connection among the `Tasks`:

```{code-cell} ipython3
flow.get_graphviz()
```

```{code-cell} ipython3
# matplotlib version based on networkx
#flow.plot_networkx(with_edge_labels=True);
```

The flow contains three independent groups of tasks, one group per each value of `ecut` specified in `params`.

```{code-cell} ipython3
for work in flow:
    for task in work:
        print(task.pos_str, "uses ecut:", task.input["ecut"])
```

```{code-cell} ipython3
flow.get_vars_dataframe("ecut")
```

Each group represents a `Workflow` and consists of one `ScfTask`(red circle) that solves the `KS` equations self-consistently producing a `WFK` file that will be used by the two children (`PhononTasks` - blue circles)
to compute the first-order change of the wavefunctions due to one of the *irreducible* atomic pertubations.

Note that `phonon_conv_flow` invokes Abinit under the hood to get the list of irreducible perturbations 
and uses this information to build the flow.
This explains why we have two `PhononTasks` per $q$-point instead of the total number of phonon modes that 
equals $3*N_{atom}=6$.

Perhaps a table with the values of the input variables associated to the DFPT perturbation will help.
`None` means that the variable is not defined in that particular input.

```{code-cell} ipython3
flow.get_vars_dataframe("rfphon", "rfatpol", "rfdir", "qpt", "kptopt")
```

If the meaning of these variables is not clear, you can consult the [Abinit documentation](https://docs.abinit.org)
e.g. [the documentation of the rfatpol input variable](https://docs.abinit.org/variables/dfpt/#rfatpol)
or access the documentation directly from python with: 

```{code-cell} ipython3
abilab.docvar("rfatpol")
```

Now we can generate the `flow_alas_ecut` directory with the input files by executing 
the `lesson_dfpt.py` script.
Then use the `abirun.py` script to launch the entire calculation with:

    abirun.py flow_alas_ecut_conv scheduler
    
You will see that all `PhononTasks` will be executed in parallel on your machine.

<div class="alert alert-warning">
Please make sure that AbiPy is properly configured by running abicheck --with flow
</div>

If you prefer to skip this part, you may want to jump to the next section, that presents the post-processing of the results.
Note that the output files are already available in the repository so it is also possible to try 
the AbiPy post-processing tools without having to run the flow.

+++

## Convergence study at $\Gamma$

There are several output files located inside the `outdata` directories:

```{code-cell} ipython3
!find flow_alas_ecut_conv/ -name "*_DDB"
```

Remember that our goal is to analyze the convergence of the phonon frequencies at $\Gamma$ 
as function of `ecut`.
So we are mainly interested in the DDB files located in the `outdata` directories 
of the `PhononWorks` (`w0/outdata`, `w1/outdata`, `w2/outdata`).
These are indeed the DDB files with all the information needed to reconstruct the 
dynamical matrix at $\Gamma$ and to compute the phonon frequencies (AbiPy calls `mrgddb`
to merge the DDB files when all the perturbations in the `PhononWork` have been computed).

The code below tells our robot that we would like to analyze all the DDB files 
located in the output directories of the works:

```{code-cell} ipython3
robot = abilab.DdbRobot.from_dir_glob("./flow_alas_ecut_conv/w*/outdata/")
robot
```

For more examples on the use of DDB and robots, see the 
[DDB notebook](https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/ddb.ipynb).
Now we ask the robot to call `anaddb` to compute the phonon frequencies at $\Gamma$ for all DDBs 
and return a pandas `DataFrame`:

```{code-cell} ipython3
data_gamma = robot.get_dataframe_at_qpoint((0, 0, 0))
```

The `DataFrame` is a dict-like object whose keys are the name of the colums in the table

```{code-cell} ipython3
print(data_gamma.keys())
```

where `mode-i` is the frequency in eV of the i-th phonon mode.

We are mainly interested in the convergence of the phonon frequencies versus `ecut` so we filter these columns with:

```{code-cell} ipython3
data_gamma = data_gamma[["ecut"] + [k for k in data_gamma if k.startswith("mode")]]
data_gamma
```

and we get some statistics about our data with:

```{code-cell} ipython3
data_gamma.describe()
```

Pandas tables are extremly powerful and the `describe` method already gives some useful info 
about the convergence of the phonon modes. 
Sometimes, however, we would like to visualize the data to have a better understanding of what's happening:

```{code-cell} ipython3
data_gamma.plot(x="ecut", y="mode3", style="-o");
```

Let's plot all the modes in different subplots with:

```{code-cell} ipython3
data_gamma.plot(x="ecut", y=[k for k in data_gamma if k.startswith("mode")], subplots=True, style="-o");
```

This convergence study at $\Gamma$ thus reveals that our pseudos require 
an `ecut` >= 6 Ha to get reasonably converged phonon frequencies at $\Gamma$.
In what follows, we assume that also the modes at the other $q$-points present a similar
convergence behaviour and we use `ecut` = 6 Ha to keep the computational cost low. 

+++

For a quick introduction to Pandas, see:

* [Pandas cookbook](https://github.com/jvns/pandas-cookbook)
* [Pandas cookbook Chapter 1 Reading data from a csv file](http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/master/cookbook/Chapter%201%20-%20Reading%20from%20a%20CSV.ipynb)

+++

## Phonon band structure of AlAs

Now we are finally ready for the calculation of the vibrational spectrum of $AlAs$.
We already managed to run DFPT calculations at $\Gamma$ with different values of `ecut` and the
steps required to get a full band structure are not that different, provided that 
the following differences are taken into account:

- we need the dynamical matrix $D(q)$ on a homogeneous mesh so that it is possible to calculate $D(R)$
  in anaddb via Fourier transform and then phonon frequencies for arbitrary q-points via Fourier interpolation
  
- $AlAs$ is a polar semiconductor so we need to include the LO-TO splitting for $q \rightarrow 0$ that, in turns,
  requires the DFPT computation of the Born effective charges and of the dielectric constant.


In AbiPy, these concepts are translated in an easy-to-use API in which you pass an initial `AbinitInput` object,
you specify the q-mesh for phonons in terms of `ph_nqpt` and activate the computation of the 
Born effective charges with the boolean flag `with_becs`.

Let's have a look at the code (as usual there are more comments than lines of code):

```{code-cell} ipython3
from lesson_dfpt import build_flow_alas_phonons
abilab.print_source(build_flow_alas_phonons)
```

We can finally construct the flow with:

```{code-cell} ipython3
flow_phbands = build_flow_alas_phonons(options=None)
```

and visualize the connections with:

```{code-cell} ipython3
flow_phbands.get_graphviz()
#flow_phbands.plot_networkx();
```

Note that there are a lot of things happening under the hood here.

First of all, AbiPy generates `PhononTasks` only for the $q$-points in the 
irreducible wedge of the Brillouin zone corresponding to `ph_ngqpt`.
Moreover, for a given $q$-point, only the irreducible atomic perturbations are explicitly computed
since the other atomic perturbations can be reconstructed by symmetry.
Fortunately you do not have to care about all these technical details as AbiPy and Abinit 
will automate the whole procedure.

Remember that the $q$-point mesh cannot be chosen arbitrarily
since all $q$ wavevectors should connect two $k$ points of the grid used for the electrons.

It is also worth stressing that the computational cost of the DFPT run depends on the q-point 
since only those symmetries that preserve the q-point as well as the direction of the perturbation 
can be employed (calculations at $\Gamma$ are therefore much faster than other q-points).

```{code-cell} ipython3
flow_phbands.get_vars_dataframe("rfphon", "rfatpol", "rfdir", "qpt", "kptopt")
```

Now we can generate the directories and the input files of the `Flow`.
Change lesson_dfpt.py so that the build_flow_alas_phonons function is called in main
instead of build_flow_alas_ecut_conv. 
Run the script to generate the directory with the flow.
Finally, use 

    abirun.py flow_alas_phonons scheduler
    
to launch the entire calculation.

+++

## Post-processing the results

Our flow is completed and we have the final DDB file with all the $q$-points and all the independent atomic perturbations. 
Let's open this DDB file with:

```{code-cell} ipython3
ddb = abilab.abiopen("flow_alas_phonons/outdata/out_DDB")
print(ddb)
```

The `DdbFile` object provides an easy-to-use interface that invokes `anaddb` to post-process
the data stored in the DDB file.

`anacompare_phdos`, for example, computes the phonon DOS with different $q$-meshes.
Each mesh is defined by a single integer, `nqsmall`, that gives the number of 
divisions used to sample the smallest reciprocal lattice vector. 
The number of divisions along the other directions are chosen so that proportions are preserved:

```{code-cell} ipython3
c = ddb.anacompare_phdos(nqsmalls=[8, 10, 12, 14, 16])
```

```{code-cell} ipython3
c.plotter.combiplot();
```

A 16x16x16 $q$-mesh with the tethraedron method gives a well converged phonon DOS.

+++

To function `anaget_phbst_and_phdos_files` allows one to compute the phonon band structure on an automatically defined $q$-path as well as the the phonon DOS:

```{code-cell} ipython3
phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(ndivsm=10, nqsmall=16, lo_to_splitting=False)

# Extract the phonon bands and the phonon DOS from phbst_file and phdos_file
phbands = phbst_file.phbands 
phdos = phdos_file.phdos
```

Let's plot the bands with matplotlib:

```{code-cell} ipython3
phbands.plot();
```

and the high-symmetry q-path:

```{code-cell} ipython3
phbands.qpoints.plot();
```

Do you see the two strange dips for the highest phonon band, at the $\Gamma$ point?
They are due to the lack of LO-TO splitting for the ANADDB treatment of the first list of vector. 
See also the discussion in the [second DFPT lesson](https://docs.abinit.org/tutorial/rf2/index.html).

+++

For years, Abinit users had to patch manually the output frequencies to include the LO-TO splitting.
These days are finally gone and we can plot the LO-TO splitting with AbiPy by just setting
lo_to_splitting=True`:

```{code-cell} ipython3
phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(ndivsm=10, nqsmall=16, lo_to_splitting=True)

# Extract the phonon bands and the phonon DOS from phbst_file and phdos_file
phbands = phbst_file.phbands 
phdos = phdos_file.phdos
```

```{code-cell} ipython3
#phbst_file.phbands.qpoints.plot(); 
```

The band structure plot now correctly shows the non-analytical behaviour around $\Gamma$:

```{code-cell} ipython3
phbands.plot();
```

<div class="alert alert-warning">
`lo_to_splitting=True` works only when the DDB contains the Born effective charges and the dielectric constant,
that must be computed in the Abinit run.
</div>

+++

To plot bands and DOS on the same figure:

```{code-cell} ipython3
phbands.plot_with_phdos(phdos);
```

The `PhdosFile` contains the phonon frequencies, the displacement vectors
as well as the decomposition of the total DOS in terms of the contributions due to 
the different types of atom in the unit cell.
This means that one can plot the type-projected phonon DOS with:

```{code-cell} ipython3
phdos_file.plot_pjdos_type(title="AlAs type-projected phonon DOS");
```

and it is even possible to plot fatbands and type-projected DOSes on the same figure with:

```{code-cell} ipython3
phbands.plot_fatbands(phdos_file=phdos_file);
```

The highest frequency modes have a strong Al-character while the low frequency modes originate from As. 
This behaviour is somehow expected. Could you explain it in terms of a simple physical model?   

+++

## Macroscopic dielectric tensor and Born effective charges

Our calculations includes the response of the system to an external electric field.
The code below extracts the macroscopic dielectric tensor (`emacro`)
and the Born effective charges (`becs`) from the DDB file:

```{code-cell} ipython3
emacro, becs = ddb.anaget_emacro_and_becs()
```

```{code-cell} ipython3
emacro
```

```{code-cell} ipython3
becs
```

As explained in the references, the Born effective charges must fulfill 
the charge neutrality sum-rule.
This rule is usually broken due to the discretization introduced by the FFT mesh, and `anaddb` will enforce it if `chneut` is set to 1 (default behaviour). Let's check it out!

```{code-cell} ipython3
print(becs)
```

Let's repeat the same calculation but without enforcing the sum-rule:

```{code-cell} ipython3
emacro, becs_chneut0 = ddb.anaget_emacro_and_becs(chneut=0)
print(becs_chneut0)
```

## Thermodynamic properties within the harmonic approximation

The thermodynamic properties of an ensemble of non-interacting phonons can be 
expresses in terms of integrals of the DOS.

\begin{equation} %\label{eq:helmholtz}
\Delta F = 3nNk_BT\int_{0}^{\omega_L}\text{ln}\left(2\text{sinh}\frac{\hbar\omega}{2k_BT}\right)g(\omega)d\omega
\end{equation}

\begin{equation} %\label{eq:free_en}
\Delta E = 3nN\frac{\hbar}{2}\int_{0}^{\omega_L}\omega\text{coth}\left(\frac{\hbar\omega}{2k_BT}\right)g(\omega)d\omega
\end{equation}

\begin{equation} %\label{eq:c_v}
C_v = 3nNk_B\int_{0}^{\omega_L}\left(\frac{\hbar\omega}{2k_BT}\right)^2\text{csch}^2\left(\frac{\hbar\omega}{2k_BT}\right)g(\omega)d\omega
\end{equation}


\begin{equation} %\label{eq:entropy}
S = 3nNk_B\int_{0}^{\omega_L}\left(\frac{\hbar\omega}{2k_BT}\text{coth}\left(\frac{\hbar\omega}{2k_BT}\right) - \text{ln}\left(2\text{sinh}\frac{\hbar\omega}{2k_BT}\right)\right)g(\omega)d\omega,
\end{equation}

where $k_B$ is the Boltzmann constant.

This should represent a reasonable approximation especially in the low temperature 
regime in which anharmonic effects can be neglected.   

Let's plot the vibrational contributions thermodynamic properties as function of $T$: 

```{code-cell} ipython3
phdos.plot_harmonic_thermo();
```

## Exercises

* Our first phonon band structure has been computed with a (4, 4, 4) $k$-mesh 
  for the electrons and a (2, 2, 2) $q$-mesh for phonons. 
  You may try to increase the density of $k$-points/$q$-points
  to see if this change affects the final results.
        
* Why do you get an error from AbiPy if you try `ngkpt` = (4, 4, 4,) and `ngqpt` = (3, 3, 3)?

+++
