#! /bin/env bash

set -o errexit
set -o nounset

: ${CMSSW_BASE:?"CMSSW_BASE is not set!  Run cmsenv!"}

PRID=$1

git pull origin master
git checkout origin/master
cd $CMSSW_BASE/src
scram b -j 4
cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test
TAG=HEAD
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG &> nano_$TAG'_data.log'
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG isMC=True &> nano_$TAG'_mc.log'

git fetch origin pull/$PRID/head:TEST_PR$PRID
git checkout TEST_PR$PRID
cd $CMSSW_BASE/src
scram b -j 4
cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test
TAG=PR$PRID
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG &> nano_$TAG'_data.log'
cmsRun run_nano_cfg.py reportEvery=10 tag=$TAG isMC=True &> nano_$TAG'_mc.log'

rm -rf validation
mkdir $TAG

python validate_nano.py testBParkNANO_data_HEAD.root testBParkNANO_data_$TAG.root
python time_analysis.py nano_HEAD_data.log --tag=HEAD
python time_analysis.py nano_$TAG'_data.log' --tag=$TAG
mv validation $TAG/validation_data

python validate_nano.py testBParkNANO_mc_HEAD.root testBParkNANO_mc_$TAG.root
python time_analysis.py nano_HEAD_mc.log --tag=HEAD
python time_analysis.py nano_$TAG'_mc.log' --tag=$TAG
mv validation $TAG/validation_mc

git checkout master
git branch -D TEST_PR$PRID
