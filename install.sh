#!/bin/bash
set -e  # exit on first error

pip install .
conda install graphviz -c conda-forge --yes
