#!/bin/sh

source ../install/activate.sh

#usage: ./run-many-cosmics-script.sh <start run number> <end run number> <events>

#do simulation in small steps
for (( c=$1; c<=$2; c++ ))
do
  run=`printf %06i $c`
  echo '/mu3e/run/setRunNumber ' $c  > run.mac
  echo '/run/beamOn ' $3 >> run.mac
  cp -f digi_cosmic.json digi.json
  cp -f digi_cosmic.json `echo 'data/mu3e_run_'$run'.json'`
  echo 'digi json copied!'
  infile=`echo 'data/mu3e_run_'$run'.root'`
  outfile=`echo 'data/mu3e_run_'$run'_trirec_cosmic.root'`
  mu3eSim run.mac
  mu3eTrirec `echo '--input '$infile' --output '$outfile' --conf trirec_cosmic.conf --cosmic'`
  rm $infile
done

