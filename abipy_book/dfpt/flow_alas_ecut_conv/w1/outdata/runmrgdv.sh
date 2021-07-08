#!/bin/bash
cd /Users/gmatteo/git_repos/abitutorials/abitutorials/dfpt/flow_alas_ecut_conv/w1/outdata
# OpenMp Environment
export OMP_NUM_THREADS=1
# Commands before execution
source ~/env.sh

mpirun -n 1 mrgdv < /Users/gmatteo/git_repos/abitutorials/abitutorials/dfpt/flow_alas_ecut_conv/w1/outdata/mrgdvdb.stdin > /Users/gmatteo/git_repos/abitutorials/abitutorials/dfpt/flow_alas_ecut_conv/w1/outdata/mrgdvdb.stdout 2> /Users/gmatteo/git_repos/abitutorials/abitutorials/dfpt/flow_alas_ecut_conv/w1/outdata/mrgdvdb.stderr
