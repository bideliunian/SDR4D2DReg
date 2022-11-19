#!/bin/bash

for args in `seq 1 100`;
do
qsub wgsir_s2d.PBS -v "args=$args"
done
