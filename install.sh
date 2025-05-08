#!/bin/bash
set -e  # exit on first error

pip install -r requirements.txt
conda install graphviz -c conda-forge --yes
