#!/bin/bash

for args in `seq 1 100`;
do
qsub wgsir_d2d.PBS -v "args=$args"
done
