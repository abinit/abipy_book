# Introduction

Welcome to the AbiPy Jupyter Book!
This book contains notebook-based documentation for [AbiPy](https://github.com/abinit/abipy).
This augments our Sphinx-based [documentation](http://pythonhosted.org/abipy/) with jupyter notebooks 
containing interactive tutorials and examples.

Additional examples are available on the:

* [AbiPy plot gallery](http://abinit.github.io/abipy/gallery/index.html)
* [AbiPy flow gallery](http://abinit.github.io/abipy/flow_gallery/index.html)
* [Matgenb website](https://matgenb.materialsvirtuallab.org)

Beginners should typically start with the three notebooks on Abipy objects,
continue with the notebook on the GSR.nc output file (and possibly some of those focusing on other output files),
and then examine the flows notebook.


+++

## AbiPy Objects

* [Structure](structure.ipynb): The crystalline structure
* [Abinit Input](abinit_input.ipynb): How to create an `AbinitInput`
* [Input Factories](input_factories.ipynb): Generating input files for high-throughput calculations

+++

## Output files supported by AbiPy

* [GSR.nc](gsr.ipynb): File with ground-state results produced by SCF/NSCF runs.  
* [HIST.nc](hist.ipynb): File produced by structural relaxations and MD runs    
* [FATBANDS.nc](efatbands.ipynb): Plot electronic fatbands and L-projected DOS
* [DDB](ddb.ipynb): Tools to analyze phonons      
* [SIGRES.nc](sigres.ipynb): How to analyze the results of $GW$ calculations 
* [MDF.nc](mdf.ipynb): How to analyze the results of Bethe-Salpeter calculations 
* [Lobster](lobster.ipynb): How to analyze the output files produced by Lobster
* [ABIWAN.nc](abiwan.ipynb): How to analyze the Wannier90 wout file and the netcdf file produced by Abinit

+++

## AbiPy Workflows

* [Flows](flows.ipynb): How to automate calculations with `Flows`, `Works` and `Tasks`

This notebook is complemented with the documentation on the [TaskManager Configuration](http://abinit.github.io/abipy/workflows/taskmanager.html): How to specify options in `manager.yml` and `scheduler.yml`  

+++

## Abinit + AbiPy Lessons

This section discusses flows and vizualisation tools for the same topics 
as some standard ABINIT tutorials.
Usually, the corresponding ABINIT tutorial has to be followed first. 
Sometimes there is also a large overlap with some of the previous AbiPy tutorials.

* [The H<sub>2</sub> molecule](base1/lesson_base1.ipynb)
* [Crystalline silicon](base3/lesson_base3.ipynb)
* [Phonons, dielectric tensor and Born effective charges from DFPT](dfpt/lesson_dfpt.ipynb)
* [$G_0W_0$ band structure](g0w0/lesson_g0w0.ipynb)
* [Bethe-Salpeter equation and excitonic effects](bse/lesson_bse.ipynb)
* [E-PH self-energy and T-dependent band structures](sigeph/lesson_sigeph.ipynb)
* [Phonon linewidths and Eliashberg function of Al](eph_al/lesson_eph.ipynb)
* [Elastic and piezoelectric properties](elastic/lesson_elastic.ipynb)

+++

## Miscellaneous

* [Useful links](links.ipynb)
