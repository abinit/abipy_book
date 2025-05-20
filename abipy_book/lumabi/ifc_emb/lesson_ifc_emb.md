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

# Phonons using the IFC embedding

This page shows how to obtain the phonon modes of a defect system using the IFC embedding approach. This technique allows one to obtain both defect-related modes and pristine-like modes on the same footing. It is first shown how to obtain the pristine and defect phonons with an AbiPy workflow using the Sr[Li$_2$Al$_2$O$_2$N$_2$]:Eu$^{2+}$ example. Then the IFC embedding is applied. 

Note that the simulation parameters are not converged, and are only provided for illustrative purposes.

## Pristine phonons

The python script below creates an AbiPy workflow that computes the phonons of the pristine system, using DFPT. It takes the primitive structure of the Sr[Li$_2$Al$_2$O$_2$N$_2$] system, and computes the phonons on a 2x2x2 q-mesh. 

```python
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.core.structure as structure
from abipy import flowtk

def make_scf_input():
    """
    This function constructs the input file for the GS calculation:
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    pseudodir='pseudos'

    pseudos = ('Eu.xml',
               'Sr.xml',
               'Al.xml',
               'Li.xml',
               'O.xml',
               'N.xml')
    stru = structure.Structure.from_file("SALON_prim.cif")

    # Initialize the input
    gs_inp = abilab.AbinitInput(stru, pseudos=pseudos,pseudo_dir=pseudodir)

    # Set variables
    gs_inp.set_vars(
        ecut=10,
        pawecutdg=20,
        diemac=5,
        nstep=100,
        tolvrs=1e-15,
        nband=60,
        nbdbuf=10,
                        )

    gs_inp.set_kmesh(ngkpt=[2,2,2],shiftk=[[0.5,0.5,0.5]])

    return gs_inp

def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")
    flow = flowtk.Flow(options.workdir, manager=options.manager)
    scf_input = make_scf_input()
    phonon_work = flowtk.works.PhononWork.from_scf_input(scf_input,qpoints=[2,2,2],is_ngqpt=True, with_becs=True,ddk_tolerance={"tolwfr":1e-10})
    flow.register_work(phonon_work)
    return flow

@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)

if __name__ == "__main__":
    sys.exit(main())
```

## Defect phonons
The script below creates an AbiPy workflow that computes the phonons of the defect system, using finite-difference. It takes the defect structure of the Sr[Li$_2$Al$_2$O$_2$] system from the $\Delta$SCF calculation, and computes the phonons at the $\Gamma$ point. 

```python
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.core.structure import Structure
from abipy import flowtk
from abipy.flowtk.abiphonopy import PhonopyWork

def make_scf_input():
    # extract the structure from the deltaSCF calculation.
    structure = Structure.from_file("flow_deltaSCF/w0/t0/outdata/out_GSR.nc")
    pseudodir='pseudos'

    pseudos = ('Eu.xml',
               'Sr.xml',
               'Al.xml',
               'Li.xml',
               'O.xml',
               'N.xml')
 
    gs_scf_inp = abilab.AbinitInput(structure, pseudos=pseudos,pseudo_dir=pseudodir)
    gs_scf_inp.set_vars(ecut=10,
                        pawecutdg=20,
                        chksymbreak=0,
                        diemac=5,
                        prtwf=-1,
                        nstep=100,
                        toldfe=1e-6,
                        chkprim=0,
                        )

    # Set DFT+U and spinat parameters according to chemical symbols.
    symb2luj = {"Eu": {"lpawu": 3, "upawu": 7, "jpawu": 0.7}}

    gs_scf_inp.set_usepawu(usepawu=1, symb2luj=symb2luj)

    n_val = gs_scf_inp.num_valence_electrons
    n_cond = round(20)

    spin_up_gs = f"\n{int((n_val - 7) / 2)}*1 7*1 {n_cond}*0"
    spin_up_ex = f"\n{int((n_val - 7) / 2)}*1 6*1 0 1 {n_cond - 1}*0"
    spin_dn = f"\n{int((n_val - 7) / 2)}*1 7*0 {n_cond}*0"

    nsppol = 2
    shiftk = [0, 0, 0]
    ngkpt = [1, 1, 1]

    # Build SCF input for the ground state configuration.
    gs_scf_inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, [spin_up_gs, spin_dn])
    # Build SCF input for the excited configuration.
    exc_scf_inp = gs_scf_inp.deepcopy()
    exc_scf_inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, [spin_up_ex, spin_dn])

    return gs_scf_inp


def build_flow(options):

    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    # Build input for GS calculation
    scf_input = make_scf_input()
    
    # We only need the phonons at Gamma point, so we set the supercell to [1,1,1].
    work = PhonopyWork.from_gs_input(scf_input, scdims=[1, 1, 1])
    flow.register_work(work)

    return flow


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())

```

## IFC Embedding

The following code cells will show how to use the two previous calculations to obtain a new Phonopy object containing the embedded phonons. 

First, we need to load the DDB file (pristine phonons) and the phonopy object of the defect system. 
```{code-cell}
from abipy.dfpt.converters import ddb_ucell_to_phonopy_supercell
from abipy.embedding.embedding_ifc import Embedded_phonons
from pymatgen.io.phonopy import get_pmg_structure
from abipy.core.kpoints import kmesh_from_mpdivs
from abipy.abilab import abiopen
import phonopy
import numpy as np


ddb_pristine=abiopen("../lineshape/flow_phonons/w0/outdata/out_DDB")
# Phonons of the unit cell bulk, computed on a 2x2x2 q-mesh (abinit DDB file)

ph_defect=phonopy.load(supercell_filename="../lineshape/flow_phonons_doped/w0/outdata/POSCAR",
                       force_sets_filename="../lineshape/flow_phonons_doped/w0/outdata/FORCE_SETS")
# Phonons obtained with defect supercell of 36 atoms (same than the Delta SCF supercell),
# obtained with finite difference with Phonopy. 
```

The folding procedure is done as follows. One provides a DDB object to the `ddb_ucell_to_phonopy_supercell` function. As the name suggests, this function takes a DDB object with a given q-mesh and returns a phonopy object with a supercell size corresponding to the DDB q-mesh, and containing the phonons modes at the Gamma point. Before calling this function, it is possible to use anaddb to obtain an interpolated DDB object, and we show that in the next code cell. 

```{code-cell}
sc_size=[2,2,4]
qpts=kmesh_from_mpdivs(mpdivs=sc_size,shifts=[0,0,0],order="unit_cell")
ddb_pristine_inter=ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts)
ph_pristine=ddb_ucell_to_phonopy_supercell(ddb_pristine_inter)
```
Explicitly, we use first a DDB file for the 18-atoms primitive structure computed 2x2x2 q-mesh, and then we use the `anaget_interpolated_ddb` function to obtain a new DDB object with an interpolated 2x2x4 q-mesh. The specific list of q-points have to be provided, and  we use the `kmesh_from_mpdivs` function to obtain this list. Finally, we use the `ddb_ucell_to_phonopy_supercell` function to obtain a phonopy object with a 2x2x4 supercell size structure (containing 288 atoms), and where the pristine IFCs were correctly folded at the Gamma point. 

We now have to combine the two phonopy objects, the one with the defect and the one with the pristine system. The `Embedded_phonons` class takes care of this. We describe below how the algorithm works.

The `Embedded_phonons.from_phonopy_instances` method combines pristine and defect phonon data to create an embedded phonon model. The process begins by converting the phonopy structures into `pymatgen` `Structure` objects and applying the necessary defect operations, such as atom substitutions, removals, or additions. Structures are cleaned and centered to make mapping between the pristine and defect supercells straightforward, ensuring that each atom in the defect supercell can be matched to its counterpart in the pristine supercell.

Once the mapping is established, the method extracts the interatomic force constants (IFCs) from both the pristine and defect phonopy objects. 

The embedding procedure is based on a spatial cutoff radius that can be set automatically or specified by the user. Within a sphere of this radius, the pristine IFCs are replaced by those from the defect calculation, ensuring that the local environment around the defect is accurately described. Outside this sphere, the pristine IFCs are retained to preserve the bulk-like character of the system. Finally, the method re-enforces the Acoustic Sum Rule (ASR).

Details on the IFCs modifications can be displayed by setting 'verbose=True'.
In the end, the method returns a new phonopy object containing the embedded phonons, which can be used for further calculations, such as lineshape simulations.

We show below a practical example of how to use the `Embedded_phonons` class.

This block of code prepares the structural informations needed for the mapping. 
In order to do that, we need to provide the defect structure without defect-induced relaxation and the coordinates of the defect in the two supercells.

```{code-cell}
# We need first to create the defect structure without relax, 
# this is the structure used to compute the defect phonons with finite difference.
structure_defect_wo_relax=ddb_pristine.structure.copy()
structure_defect_wo_relax.make_supercell([1,1,2])
structure_defect_wo_relax.replace(0,'Eu')

# index of the sub. = 0 (in defect structure), this is found manually
idefect_defect_stru=0
main_defect_coords_in_defect=structure_defect_wo_relax.cart_coords[idefect_defect_stru]

# index of the sub. = 0 (in pristine structure), this is found manually
from pymatgen.io.phonopy import get_pmg_structure
idefect_pristine_stru=0
main_defect_coords_in_pristine=get_pmg_structure(ph_pristine.supercell).cart_coords[idefect_pristine_stru]
```

We now call the embedding algorithm, which creates a phonopy object containing the embedded phonons. 

```{code-cell}
help(Embedded_phonons.from_phonopy_instances)
```

```{code-cell}
emb_ph=Embedded_phonons.from_phonopy_instances(phonopy_pristine=ph_pristine,
                                               phonopy_defect=ph_defect,
                                               structure_defect_wo_relax=structure_defect_wo_relax,
                                               main_defect_coords_in_defect=main_defect_coords_in_defect,
                                               main_defect_coords_in_pristine=main_defect_coords_in_pristine,
                                               substitutions_list=[[idefect_pristine_stru,"Eu"]],
                                               cut_off_mode="auto",
                                       ) 
```


If the structure mapping fails (like in the code below), the code will tell the user that the mapping did not work:  
```{code-cell}
emb_ph_failed=Embedded_phonons.from_phonopy_instances(phonopy_pristine=ph_pristine,
                                               phonopy_defect=ph_defect,
                                               structure_defect_wo_relax=structure_defect_wo_relax,
                                               main_defect_coords_in_defect=main_defect_coords_in_defect,
                                               main_defect_coords_in_pristine=main_defect_coords_in_pristine+np.array([0.1,0.0,0.0]),
                                               substitutions_list=[[idefect_pristine_stru,"Eu"]],
                                               cut_off_mode="auto",verbose=False
                                       ) 
```



