#!/bin/bash

for args in `seq 1 99`;
do
qsub wgsir_d2v.PBS -v "args=$args"
done
