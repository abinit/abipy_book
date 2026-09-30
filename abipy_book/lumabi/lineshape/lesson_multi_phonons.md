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

# Lineshape (multi-phonon)

This tutorial shows how to obtain the luminescence lineshape within a multi-phonon mode model.

```{code-cell}
import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
import numpy as np
import phonopy

from pymatgen.io.phonopy import get_pmg_structure
from abipy.abilab import abiopen
from abipy.dfpt.converters import ddb_ucell_to_phonopy_supercell
from abipy.embedding.embedding_ifc import Embedded_phonons
from abipy.lumi.deltaSCF import DeltaSCF
from abipy.core.kpoints import kmesh_from_mpdivs
from abipy.lumi.lineshape import Lineshape
```

We first load the ab-initio phonon results.
The inputs used to generate these calculations are shown in the next [section](../ifc_emb/lesson_ifc_emb.md).

```{code-cell}
ddb_pristine = abiopen("../workflows_data/flow_phonons/w0/outdata/out_DDB")

# Phonons of the unit cell bulk, computed on a q-mesh (abinit DDB file)

ph_defect = phonopy.load(supercell_filename="../workflows_data/flow_phonons_doped/w0/outdata/POSCAR",
                         force_sets_filename="../workflows_data/flow_phonons_doped/w0/outdata/FORCE_SETS")

# Phonons obtained with defect supercell of 36 atoms (same than the Delta SCF supercell),
# obtained with finite difference with Phonopy.
```

Then we create a `DeltaSCF` object from the results of a LumiWork workflow:

```{code-cell}
files = [
    "../workflows_data/flow_deltaSCF/w0/t2/outdata/out_GSR.nc",
    "../workflows_data/flow_deltaSCF/w0/t3/outdata/out_GSR.nc",
    "../workflows_data/flow_deltaSCF/w0/t4/outdata/out_GSR.nc",
    "../workflows_data/flow_deltaSCF/w0/t5/outdata/out_GSR.nc",
]

results = DeltaSCF.from_four_points_file(files)
```

## Simplest case, no IFCs embedding

We create a lineshape object by combining the $\Delta$SCF and phonon calculations.
In this example, the coordinates of the defect are the same in the two supercells,
but this might not always be the case,
which is why these values must be provided explicitly.
Notice that, in this first simple example, the phonon supercell has the same size as the $\Delta$SCF supercell,
so either the displacements or the forces can be used.

```{code-cell}
# get_pmg_structure(ph_defect.supercell) # to inspect the coords of the defect in phonon scell

lineshape = Lineshape.from_phonopy_phonons(
    E_zpl=results.E_zpl(),
    phonopy_ph=ph_defect,
    dSCF_structure=results.structure_gs(),
    use_forces=False, # we choose the displacements.
    dSCF_displacements=results.diff_pos(),
    dSCF_forces=results.forces_gs,
    coords_defect_dSCF=np.array([3.9795, 0.0000, 1.5920]),
    coords_defect_phonons=np.array([3.9795, 0.0000, 1.5920])
)

help(Lineshape.from_phonopy_phonons)
```

We can now plot the Huang-Rhys spectral function $S(\hbar\omega)=\sum_{\nu}S_{\nu}\delta(\hbar\omega-\hbar\omega_{\nu})$,
optionally with colors showing the degree of localization of each mode:

```{code-cell}
lineshape.plot_spectral_function(with_S_nu=True);
```

```{code-cell}
lineshape.plot_spectral_function(with_local_ratio=True);
```

The total Huang-Rhys factor computed with the multi-phonon mode model can be compared with the one obtained with the one-dimensional model:

```{code-cell}
print(f"multi phonon Huang-Rhys factor  = {np.round(lineshape.S_tot(),3)}")
print(f"one  phonon Huang-Rhys factor  = {np.round(results.S_em(),3)}")
```

The final luminescence emission spectrum, computed with the generating function approach, can be plotted at any temperature with:

```{code-cell}
x, y = lineshape.L_hw(T=300, w=5) # w is the width of the gaussian used to smooth the spectrum.
fig, ax = plt.subplots(figsize=(6,3))
ax.plot(x, y)
ax.set_xlabel("Energy (eV)")
ax.set_ylabel("Emission Intensity (a.u.)")
ax.set_xlim(1.2,1.8)
```

## With IFCs embedding

We first need to create a phonopy object with embedded IFCs
(see the next tutorial for more details).

The first block of code interpolates the DDB file on the desired $q$-mesh (here 2x2x4).
It then folds these phonons from the unit cell on a 2x2x4 $q$-mesh to the phonons
of the corresponding 2x2x4 supercell at the $\Gamma$ point, and saves them in phonopy format.

```{code-cell}
sc_size = [2,2,4]
qpts = kmesh_from_mpdivs(mpdivs=sc_size,shifts=[0,0,0],order="unit_cell")
ddb_pristine_inter = ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts)
ph_pristine = ddb_ucell_to_phonopy_supercell(ddb_pristine_inter)
```

This block of code prepares the structural information needed for the mapping.

```{code-cell}
# We need first to create the defect structure without relax

structure_defect_wo_relax = ddb_pristine.structure.copy()
structure_defect_wo_relax.make_supercell([1, 1, 2])
structure_defect_wo_relax.replace(0, 'Eu')

# index of the sub. = 0 (in defect structure), this is found manually
idefect_defect_stru = 0
main_defect_coords_in_defect = structure_defect_wo_relax.cart_coords[idefect_defect_stru]

# index of the sub. = 0 (in pristine structure), this is found manually
from pymatgen.io.phonopy import get_pmg_structure
idefect_pristine_stru=0
main_defect_coords_in_pristine = get_pmg_structure(ph_pristine.supercell).cart_coords[idefect_pristine_stru]
```

We now call the embedding algorithm, which creates a phonopy object containing the embedded phonons:

```{code-cell}
emb_ph = Embedded_phonons.from_phonopy_instances(
    phonopy_pristine=ph_pristine,
    phonopy_defect=ph_defect,
    structure_defect_wo_relax=structure_defect_wo_relax,
    main_defect_coords_in_defect=main_defect_coords_in_defect,
    main_defect_coords_in_pristine=main_defect_coords_in_pristine,
    substitutions_list=[[idefect_pristine_stru,"Eu"]],
    cut_off_mode="auto",verbose=False
)
```

The lineshape object can now be created with these new phonons:

```{code-cell}
lineshape_emb = Lineshape.from_phonopy_phonons(
    E_zpl=results.E_zpl(),
    phonopy_ph=emb_ph,
    dSCF_structure=results.structure_gs(),
    use_forces=True,
    dSCF_displacements=results.diff_pos(),
    dSCF_forces=results.forces_gs,
    coords_defect_dSCF=np.array([3.9795, 0.0000, 1.5920]),
    coords_defect_phonons=np.array([0,0,0]), # note that the defect coords changed!
)
```

```{code-cell}
lineshape_emb.plot_spectral_function(with_local_ratio=True);
```

```{code-cell}
x, y = lineshape_emb.L_hw(T=300, w=5) # w is the width of the gaussian used to smooth the spectrum.
fig, ax = plt.subplots(figsize=(6,3))
ax.plot(x, y)
ax.set_xlabel("Energy (eV)")
ax.set_ylabel("Emission Intensity (a.u.)")
ax.set_xlim(1.2, 1.8)
```

Notice how the spectrum has broadened with the increase of the supercell size.
This is because the computed Huang-Rhys factor is larger than the one obtained with the smaller supercell,
which illustrates the importance of converging the supercell size.
The fact that the phonon peaks are less resolved is also due to the increased number of phonon modes.
You can also play with the smoothing parameters `w` (Gaussian) and `lamb` (Lorentzian) to see how they affect the spectrum.

```{code-cell}
print(f"multi phonon Huang-Rhys factor with larger supercell and embedding = {np.round(lineshape_emb.S_tot(),3)}")
```

## Convergence with respect to the supercell size

The following block of code illustrates how to perform a convergence study with respect to the supercell size.
The supercell sizes are defined in the `sc_sizes` list.
The code loops over these sizes and computes the lineshape for each of them.
We then plot the Huang-Rhys spectral function for each supercell size.

```{code-cell}
sc_sizes=[ [1, 1, 2],
           [2, 2, 2],
           [2, 2, 4],]

lineshape_emb_list = []

structure_defect_wo_relax = ddb_pristine.structure.copy()
structure_defect_wo_relax.make_supercell([1, 1, 2])
structure_defect_wo_relax.replace(0, 'Eu')

# index of the sub. = 0 (in defect structure), this is found manually
idefect_defect_stru = 0
main_defect_coords_in_defect = structure_defect_wo_relax.cart_coords[idefect_defect_stru]

for sc_size in sc_sizes:
    print(f"supercell size: {sc_size}")
    qpts = kmesh_from_mpdivs(mpdivs=sc_size,shifts=[0,0,0],order="unit_cell")
    ddb_pristine_inter = ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts)
    ph_pristine = ddb_ucell_to_phonopy_supercell(ddb_pristine_inter)

    # index of the sub. = 0 (in pristine structure), this is found manually
    from pymatgen.io.phonopy import get_pmg_structure
    idefect_pristine_stru = 0
    main_defect_coords_in_pristine = get_pmg_structure(ph_pristine.supercell).cart_coords[idefect_pristine_stru]

    emb_ph = Embedded_phonons.from_phonopy_instances(
        phonopy_pristine=ph_pristine,
        phonopy_defect=ph_defect,
        structure_defect_wo_relax=structure_defect_wo_relax,
        main_defect_coords_in_defect=main_defect_coords_in_defect,
        main_defect_coords_in_pristine=main_defect_coords_in_pristine,
        substitutions_list=[[idefect_pristine_stru,"Eu"]],
        cut_off_mode="auto",verbose=False
    )

    lineshape_emb = Lineshape.from_phonopy_phonons(
        E_zpl=results.E_zpl(),
        phonopy_ph=emb_ph,
        dSCF_structure=results.structure_gs(),
        use_forces=True,
        dSCF_displacements=results.diff_pos(),
        dSCF_forces=results.forces_gs,
        coords_defect_dSCF=np.array([3.9795, 0.0000, 1.5920]),
        coords_defect_phonons=np.array([0,0,0])) # note that the defect coords changed!

    lineshape_emb_list.append(lineshape_emb)
```

```{code-cell}
  fig, axs = plt.subplots(figsize=(8,4))

  for i,lineshape in enumerate(lineshape_emb_list):
      x, y = lineshape.S_hbarOmega(broadening=3)
      legend = f"supercell size: {sc_sizes[i]}, S={np.round(lineshape.S_tot(),2)}"
      axs.plot(x, y, label=legend)

  axs.legend()
  axs.set_xlabel("Phonon energy (eV)")
  axs.set_ylabel(r"$S(\hbar\omega$)")
```

```{note}
For additional examples of the use of this module, see the tests in `abipy/lumi/tests/test_lineshape.py`.
```
