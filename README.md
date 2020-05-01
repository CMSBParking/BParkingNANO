# nanoAOD producer customized for BParking analysis (focus on RK/K*/phi)

## Getting started

```shell
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

[]: # (Old version: Add the latest code and model (2019Aug07) for the electron ID)
## Add the latest code and model (2019Aug07) for the electron ID 
#
#```shell
#git cms-addpkg RecoEgamma/EgammaElectronProducers
#git cms-merge-topic CMSBParking:from-CMSSW_10_2_15_2019Aug07
#git cms-addpkg RecoEgamma/ElectronIdentification
#scram b
#
# Check $CMSSW_BASE/external exists before this step (e.g. run 'scram b' to create it)
#git clone --single-branch --branch 102X_LowPtElectrons_2019Aug07 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
#
# The following step is required if running on CRAB
#mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data 


## Add the latest code and model (2019Aug07 in mvaId and Feb24-depth10 in mvaIdExtra) for the electron ID

```shell
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-addpkg RecoEgamma/ElectronIdentification
git remote add crovelli git@github.com:crovelli/cmssw.git
git cms-merge-topic crovelli:from-CMSSW_10_2_15_2020Feb24-depth10


## Add the modification needed to use post-fit quantities for electrons  

```shell
git cms-addpkg TrackingTools/TransientTrack
git cms-merge-topic -u CMSBParking:GsfTransientTracks
```

## Add the modification needed to use the KinematicParticleVertexFitter  

```shell
git cms-merge-topic -u CMSBParking:fixKinParticleVtxFitter
```

## Add the BParkingNano package and build everything

```shell
git clone git@github.com:CMSBParking/BParkingNANO.git  ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```

## To run on a test file

```shell
cd PhysicsTools/BParkingNano/test/
cmsenv 
cmsRun run_nano_cfg.py
```

## Contributing

We use the _fork and pull_ model:

fork this repository https://github.com/CMSBParking/BParkingNANO (top right _Fork button)

If you haven't done so yet, clone this repository:

```shell
git clone git@github.com:CMSBParking/BParkingNANO.git  ./PhysicsTools
```

Add your fork of the repository as remote:

```shell
git remote add mine git@github.com:`git config user.github`/BParkingNANO.git
git checkout -b ${USER}_feature_branch origin/master
```

Work on your feature, `add`, `commit`, etc. and push to your own fork

when adding a sequence or table producer, please include it in the _python/nanoBPark_cff.py_
and make sure it runs properly checking the output result (_test_BParkSequence_10215.py_ to give it a try)

```shell
git push mine feature_branch
```

Make a pull request on github
