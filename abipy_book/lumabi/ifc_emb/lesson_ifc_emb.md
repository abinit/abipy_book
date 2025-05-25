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

# IFCs Embedding

This section explains how to obtain the phonon modes of a defect system using the IFC (Interatomic Force Constant) embedding approach.
This method enables the calculation of both defect-localized and bulk-like phonon modes within a unified framework.
We illustrate the process using the Sr[Li$_2$Al$_2$O$_2$N$_2$]:Eu$^{2+}$ example, first computing pristine and defect phonons,
then performing the embedding.

```{note}
The simulation parameters below are for demonstration only and are not converged.
```

## 1. Pristine Phonons

The following script creates an AbiPy workflow to compute the phonons of the pristine system using DFPT.
It uses the primitive structure of Sr[Li$_2$Al$_2$O$_2$N$_2$] and computes phonons on a 2x2x2 q-mesh.

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
    pseudodir = 'pseudos'

    pseudos = ('Eu.xml',
               'Sr.xml',
               'Al.xml',
               'Li.xml',
               'O.xml',
               'N.xml')
    stru = structure.Structure.from_file("SALON_prim.cif")

    # Initialize the input
    gs_inp = abilab.AbinitInput(stru, pseudos=pseudos, pseudo_dir=pseudodir)

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

    gs_inp.set_kmesh(ngkpt=[2,2,2], shiftk=[[0.5,0.5,0.5]])

    return gs_inp


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(options.workdir, manager=options.manager)
    scf_input = make_scf_input()
    phonon_work = flowtk.works.PhononWork.from_scf_input(scf_input, qpoints=[2,2,2], is_ngqpt=True,
                                                         with_becs=True, ddk_tolerance={"tolwfr":1e-10})
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

## 2. Defect Phonons

This script computes the phonons of the defect system using finite differences.
It uses the supercell defect structure from the $\Delta$SCF calculation and then uses a phonopy workflow.
The supercell size is set to [1,1,1], which is equivalent to a $\Gamma$ point calculation.

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
    pseudodir = 'pseudos'

    pseudos = ('Eu.xml',
               'Sr.xml',
               'Al.xml',
               'Li.xml',
               'O.xml',
               'N.xml')

    gs_scf_inp = abilab.AbinitInput(structure, pseudos=pseudos, pseudo_dir=pseudodir)
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

## 3. IFC Embedding: Step-by-Step

The IFC embedding procedure combines the results of the two previous calculations to produce a new `Phonopy` object
that contains the phonon modes of the defect system in a large supercell.
The embedded phonons are typically created with the `Embedded_phonons.from_phonopy_instances` method:

```{code-cell}
from abipy.embedding.embedding_ifc import Embedded_phonons
help(Embedded_phonons.from_phonopy_instances)
```

### 3.1. Load Pristine and Defect Phonon Data

First, load the DDB file (pristine phonons) and the phonopy object for the defect system:

```{code-cell}
from abipy.dfpt.converters import ddb_ucell_to_phonopy_supercell
from abipy.embedding.embedding_ifc import Embedded_phonons
from pymatgen.io.phonopy import get_pmg_structure
from abipy.core.kpoints import kmesh_from_mpdivs
from abipy.abilab import abiopen
import phonopy
import numpy as np

# Load pristine phonons (DFPT, 2x2x2 q-mesh)
ddb_pristine = abiopen("../workflows_data/flow_phonons/w0/outdata/out_DDB")

# Load defect phonons (finite difference, supercell)
ph_defect = phonopy.load(
    supercell_filename="../workflows_data/flow_phonons_doped/w0/outdata/POSCAR",
    force_sets_filename="../workflows_data/flow_phonons_doped/w0/outdata/FORCE_SETS"
)
```

### 3.2. Fold Pristine IFCs to Supercell

The pristine DDB file is first interpolated using `anaget_interpolated_ddb` and then the folding procedure
is done with `ddb_ucell_to_phonopy_supercell`:

```{code-cell}
# Define the large phonon supercell size for the embedding
sc_size = [2, 2, 4]

# Generate the list of q-points corresponding to the supercell
qpts = kmesh_from_mpdivs(mpdivs=sc_size, shifts=[0, 0, 0], order="unit_cell")

# Interpolate the pristine DDB to obtain force constants on the q-point mesh
ddb_pristine_inter = ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts)

# Fold the phonons, and convert to Phonopy object.
ph_pristine = ddb_ucell_to_phonopy_supercell(ddb_pristine_inter)
```

### 3.3. Prepare Structures and Defect Mapping

To embed the IFCs, the algorithm needs to know which atom(s) are substituted or modified.
This requires:

- The defect structure *without* relaxation (i.e., before local relaxation around the defect)
- The coordinates of the defect atom in both the defect and pristine supercells

```{code-cell}
# Create the defect structure without relaxation
structure_defect_wo_relax = ddb_pristine.structure.copy()
structure_defect_wo_relax.make_supercell([1, 1, 2])
structure_defect_wo_relax.replace(0, 'Eu')

# Index and coordinates of the defect atom (manually identified)
idefect_defect_stru = 0
main_defect_coords_in_defect = structure_defect_wo_relax.cart_coords[idefect_defect_stru]

idefect_pristine_stru = 0 # (manually identified)
main_defect_coords_in_pristine = get_pmg_structure(ph_pristine.supercell).cart_coords[idefect_pristine_stru]
```

### 3.4. Run the Embedding Algorithm

The `Embedded_phonons.from_phonopy_instances` method combines the pristine and defect phonon data.
The algorithm:

- Maps atoms between the pristine and defect supercells (crucial step)
- Extracts IFCs from both systems
- Replaces pristine IFCs with defect IFCs within a cutoff radius around the defect
- Retains pristine IFCs elsewhere
- Enforces the Acoustic Sum Rule (ASR)

You can print details of the modifications by setting `verbose=True`.

```{code-cell}
emb_ph = Embedded_phonons.from_phonopy_instances(
    phonopy_pristine=ph_pristine,
    phonopy_defect=ph_defect,
    structure_defect_wo_relax=structure_defect_wo_relax,
    main_defect_coords_in_defect=main_defect_coords_in_defect,
    main_defect_coords_in_pristine=main_defect_coords_in_pristine,
    substitutions_list=[[idefect_pristine_stru, "Eu"]],
    cut_off_mode="auto"
)
```

If the structure mapping fails (e.g., due to mismatched coordinates), the code will notify you:

```{code-cell}
emb_ph_failed = Embedded_phonons.from_phonopy_instances(
    phonopy_pristine=ph_pristine,
    phonopy_defect=ph_defect,
    structure_defect_wo_relax=structure_defect_wo_relax,
    main_defect_coords_in_defect=main_defect_coords_in_defect,
    main_defect_coords_in_pristine=main_defect_coords_in_pristine + np.array([0.1, 0.0, 0.0]),
    substitutions_list=[[idefect_pristine_stru, "Eu"]],
    cut_off_mode="auto",
    verbose=False
)
```

Finally, note that the `Embedded_phonons` class provides methods for further analysis and data export.
You can use `get_gamma_freq_with_vec_abipy_fmt()` to compute the $\Gamma$-point phonon frequencies and eigenvectors in AbiPy format,
or `to_ddb()` to convert the embedded phonons back to an Abinit DDB file.
Also, you will find in `abipy.embedding.utils_ifc` function to compute the localization ratio of the phonon modes
as well as helper function to draw phonon eigenvectors with VESTA.


```{note}
For additional examples of the use of this module, especially for the embedding of different defect types (vacancy, interstitial),
you can check the tests examples located in `abipy/embedding/tests/test_embedding_ifc.py`.
```
