#!/bin/bash

KARIPATH="../3rdparty/karimaki/"
SRC="recacc_main.cpp reconstruction_accuracy.cpp ${KARIPATH}karimaki_hit.c ${KARIPATH}fitszw.c ${KARIPATH}cirpar.c ${KARIPATH}cirparw.c ${KARIPATH}wrappers.c"

ROOT=" `root-config --libs --cflags`"

echo
echo "Compiling reconstruction accuracy ..."
echo "sourcefiles: $SRC"

g++ -o Mu3eRecAcc $SRC -lf2c -DMAIN $ROOT