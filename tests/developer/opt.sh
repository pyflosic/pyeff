#!/bin/bash

# path to the pyeff src code 
export PYTHONPATH=../../../../src:$PYTHONPATH
python pyeff_optimizer_pbc.py > pyeff.log 
