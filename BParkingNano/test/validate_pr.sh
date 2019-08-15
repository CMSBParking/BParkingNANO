#! /bin/env bash

PRID=$1
if [ $# -eq 2 ]; then
    remote=$2
else
    remote=origin
fi

set -o errexit
set -o nounset

: ${CMSSW_BASE:?"CMSSW_BASE is not set!  Run cmsenv!"}

echo "Getting the latest HEAD of the common repository"
git fetch $remote master
git checkout $remote/master -b official_current_master
cd $CMSSW_BASE/src
echo "Compiling..."
scram b -j 4
echo "...Done!"
cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test
TAG=HEAD
echo "Getting reference for data..."
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG &> nano_$TAG'_data.log'
echo "Done! Now for MC..."
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG isMC=True &> nano_$TAG'_mc.log'

echo "Now merging the changes for PR #"$PRID
git fetch $remote pull/$PRID/head:TEMP_PR$PRID
git checkout official_current_master -b TEST_PR$PRID
git merge --no-edit TEMP_PR$PRID
cd $CMSSW_BASE/src
echo "Compiling..."
scram b -j 4 > PhysicsTools/BParkingNano/test/compilation_PR$PRID.log
echo "...Done!"
cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test
TAG=PR$PRID
echo "Testing on data..."
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG &> nano_$TAG'_data.log'
echo "... Done! And now on MC..."
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG isMC=True &> nano_$TAG'_mc.log'
echo "...Done! Making validation plots"

rm -rf validation
mkdir $TAG
mv compilation_PR$PRID.log $TAG/.

python validate_nano.py testBParkNANO_data_HEAD.root testBParkNANO_data_$TAG.root --plot-only-failing
python time_analysis.py nano_HEAD_data.log --tag=HEAD
python time_analysis.py nano_$TAG'_data.log' --tag=$TAG
mv validation $TAG/validation_data

python validate_nano.py testBParkNANO_mc_HEAD.root testBParkNANO_mc_$TAG.root --plot-only-failing
python time_analysis.py nano_HEAD_mc.log --tag=HEAD
python time_analysis.py nano_$TAG'_mc.log' --tag=$TAG
mv validation $TAG/validation_mc

echo "Getting rid of unused branches"
git checkout master
git branch -D TEST_PR$PRID
git branch -D TEMP_PR$PRID
git branch -D official_current_master

echo "EVERYTHING DONE! you can find the validation plots and html-based text in:" $CMSSW_BASE/src/PhysicsTools/BParkingNano/test/$TAG
