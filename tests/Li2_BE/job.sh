#!/bin/bash

# path to the pyeff src code 
export PYTHONPATH=../../src:$PYTHONPATH
# ignore python warnings 
export PYTHONWARNINGS="ignore"
python run.py > pyeff.log 
