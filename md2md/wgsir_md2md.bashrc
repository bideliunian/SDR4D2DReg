#!/bin/bash

for args in `seq 1 100`;
do
qsub wgsir_md2md.PBS -v "args=$args"
done
