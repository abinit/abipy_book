import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.core.structure import Structure
from abipy import flowtk
from abipy.flowtk.abiphonopy import PhonopyWork


def make_scf_input():
    """
    This function constructs the input file for the GS calculation:
    """
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
#    symb2spinat = {"Eu": [0, 0, 7]}
    symb2luj = {"Eu": {"lpawu": 3, "upawu": 7, "jpawu": 0.7}}

    gs_scf_inp.set_usepawu(usepawu=1, symb2luj=symb2luj)
#    gs_scf_inp.set_spinat_from_symbols(symb2spinat, default=(0, 0, 0))

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
    """
    Create a `Flow` for phonon calculations. The flow has two works.

    The first work contains a single GS task that produces the WFK file used in DFPT
    Then we have multiple Works that are generated automatically
    in order to compute the dynamical matrix on a [4, 4, 4] mesh.
    Symmetries are taken into account: only q-points in the IBZ are generated and
    for each q-point only the independent atomic perturbations are computed.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    # Build input for GS calculation
    scf_input = make_scf_input()

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
