#!/bin/bash

# path to the pyeff src code 
export ASE_EFF_PBC_COMMAND="/home/schwalbe/Programms/eff_distro_devil/eff/eff eff.cfg"
export PYTHONPATH=../../../../../src:$PYTHONPATH
python eff_optimizer_pbc.py > eff.log 
