#!/bin/bash
cd /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_phonons/w0/outdata
# OpenMp Environment
export OMP_NUM_THREADS=1
# Commands before execution
export PATH=/Users/jbouquiaux/git/abinit/build/src/98_main:$PATH

mpirun  -n 1 mrgdv < /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_phonons/w0/outdata/mrgdvdb.stdin > /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_phonons/w0/outdata/mrgdvdb.stdout 2> /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_phonons/w0/outdata/mrgdvdb.stderr
