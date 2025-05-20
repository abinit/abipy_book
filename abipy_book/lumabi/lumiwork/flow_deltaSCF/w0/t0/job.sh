#!/bin/bash
cd /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_deltaSCF/w0/t0
# OpenMp Environment
export OMP_NUM_THREADS=1
# Commands before execution
export PATH=/Users/jbouquiaux/git/abinit/build/src/98_main:$PATH

mpirun  -n 2 abinit /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_deltaSCF/w0/t0/run.abi  > /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_deltaSCF/w0/t0/run.log 2> /Users/jbouquiaux/git/lumi_book/lumi_book/lumi/flow_deltaSCF/w0/t0/run.err
