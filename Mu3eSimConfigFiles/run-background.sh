source ../install/activate.sh
echo '/mu3e/run/setRunNumber ' $1  > run.mac
echo '/run/beamOn ' $2 >> run.mac
run=`printf %06i $1`

cp -f digi_bg.json digi.json
cp -f digi_bg.json `echo 'data/mu3e_run_'$run'.json'`

mu3eSim run.mac
