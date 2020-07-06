#!/bin/sh

echo '/mu3e/run/setRunNumber ' $1  > run.mac
echo '/run/beamOn ' $2 >> run.mac
run=`printf %06i $1`
infile=`echo 'data/mu3e_run_'$run'.root'`
outfile=`echo 'data/mu3e_run_'$run'_trirec_cosmic.root'`
mu3eSim run.mac
mu3eTrirec `echo '--input ' $infile' --output '$outfile' --conf trirec_cosmic.conf --cosmic'