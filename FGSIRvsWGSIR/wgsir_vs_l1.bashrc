#!/bin/bash

for args in `seq 1 100`;
do
qsub wgsir_vs_l1.PBS -v "args=$args"
done
