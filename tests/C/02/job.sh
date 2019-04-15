#!/bin/bash

# path to the pyeff src code 
export PYTHONPATH=../../../src:$PYTHONPATH
python run_nrlmol.py > pyeff.log 
