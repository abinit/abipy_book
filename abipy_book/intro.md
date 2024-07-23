# Introduction

Welcome to the AbiPy Jupyter Book!
This book contains notebook-based documentation for [AbiPy](https://github.com/abinit/abipy).
This augments our Sphinx-based [documentation](https://abinit.github.io/abipy>) with jupyter notebooks
containing interactive tutorials and examples.

Additional examples are available on the:

* [AbiPy plot gallery](https://abinit.github.io/abipy/gallery/index.html)
* [AbiPy flow gallery](https://abinit.github.io/abipy/flow_gallery/index.html)
* [Matgenb website](https://matgenb.materialsvirtuallab.org)

Beginners should typically start with the three notebooks on Abipy objects,
continue with the notebook on the GSR.nc output file (and possibly some of those focusing on other output files),
and then examine the flows notebook.

+++

## AbiPy Objects

* [Structure](structure): The crystalline structure
* [Abinit Input](abinit_input): How to create an `AbinitInput`
* [Input Factories](input_factories): Generating input files for high-throughput calculations

+++

## Output files supported by AbiPy

* [GSR.nc](gsr): File with ground-state results produced by SCF/NSCF runs.
* [HIST.nc](hist): File produced by structural relaxations and MD runs
* [FATBANDS.nc](efatbands): Plot electronic fatbands and L-projected DOS
* [DDB](ddb): Tools to analyze phonons
* [SIGRES.nc](sigres): How to analyze the results of $GW$ calculations
* [MDF.nc](mdf): How to analyze the results of Bethe-Salpeter calculations
* [Lobster](lobster): How to analyze the output files produced by Lobster
* [ABIWAN.nc](abiwan): How to analyze the Wannier90 wout file and the netcdf file produced by Abinit

+++

## AbiPy Workflows

* [Flows](/flows.html): How to automate calculations with `Flows`, `Works` and `Tasks`

This notebook is complemented with the documentation on the
[TaskManager Configuration](https://abinit.github.io/abipy/workflows/taskmanager.html):
How to specify options in `manager.yml` and `scheduler.yml`

+++

## Abinit + AbiPy Lessons

This section discusses flows and vizualisation tools for the same topics
as some standard ABINIT tutorials.
Usually, the corresponding ABINIT tutorial has to be followed first.
Sometimes there is also a large overlap with some of the previous AbiPy tutorials.

* [The H<sub>2</sub> molecule](base1/lesson_base1)
* [Crystalline silicon](base3/lesson_base3)
* [Phonons, dielectric tensor and Born effective charges from DFPT](dfpt/lesson_dfpt)
* [$G_0W_0$ band structure](g0w0/lesson_g0w0)
* [Bethe-Salpeter equation and excitonic effects](bse/lesson_bse)
% * [E-PH self-energy and T-dependent band structures](eph_zpr/lesson_eph_zpr)
% * [Phonon linewidths and Eliashberg function of Al](eph_isotc/lesson_eph_isotc)
* [Elastic and piezoelectric properties](elastic/lesson_elastic)

+++

