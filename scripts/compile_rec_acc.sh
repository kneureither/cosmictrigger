#!/bin/bash

SRC=" karimaki_main.c karimaki_hit.c fitszw.c cirpar.c cirparw.c wrappers.c"

ROOT=" `root-config --libs --cflags`"

echo
echo "Compiling reconstruction accuracy ..."
echo "sourcefiles: $SRC"

g++ -o Mu3eRecAcc recacc_main.c reconstruction_accuracy.c $ROOT