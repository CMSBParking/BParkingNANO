## nanoAOD producer customized for BParking analysis (focus on RK/K*/phi)


## Getting started

```shell
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```
## Add the latest code and model for the electron ID 
```
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-merge-topic CMSBParking:from-CMSSW_10_2_15_LowPtElectronsID
git cms-addpkg RecoEgamma/ElectronIdentification
git clone --single-branch --branch 102X_LowPtElectrons_2019Jun28 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data # this is required if running on CRAB
```

## Add the modification needed to use post-fit quantities for electrons  
```
git cms-addpkg TrackingTools/TransientTrack
git cms-merge-topic -u CMSBParking:GsfTransientTracks
```

## Add the BParkingNano package and build everything

```
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
