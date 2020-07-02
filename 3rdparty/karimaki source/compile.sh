#!/bin/bash

SRC=" karimaki_main.c karimaki_hit.c fitszw.c cirpar.c cirparw.c wrappers.c"
ROOT=" `root-config --libs --cflags`"

echo
echo "Compiling karimaki fit..."
echo "sourcefiles: $SRC"

g++ -o karifit $SRC -lf2c -DMAIN $ROOT