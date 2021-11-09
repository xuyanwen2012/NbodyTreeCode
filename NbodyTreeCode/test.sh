#!/bin/bash

rm tmp_new.txt
rm mlpStats memStats epochStats
make
./decades_base/decades_base 1
../../MosaicSim_private/tools/mosaicrun -cc yanwen_inorder . > tmp_new.txt
cat tmp_new.txt
python3 show_mem.py



