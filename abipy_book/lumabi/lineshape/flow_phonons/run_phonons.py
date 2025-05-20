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
    gs_inp.set_vars(ecut=10,
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
