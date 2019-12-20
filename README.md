## nanoAOD producer customized for BParking analysis (focus on RK/K*/phi)


## Getting started

```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init

```
## Add the latest code and 2019Aug07 model for the electron ID 
```
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-merge-topic CMSBParking:from-CMSSW_10_2_15_2019Aug07
git cms-addpkg RecoEgamma/ElectronIdentification
scram b

# Check $CMSSW_BASE/external exists before this step (e.g. just run scram b to create it)
git clone --single-branch --branch 102X_LowPtElectrons_2019Aug07 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data

# The following step is required if running on CRAB
mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data 
```

## Add the BParkingNano package and build everything

```
git clone git@github.com:CMSBParking/BParkingNANO.git ./PhysicsTools
scram b
```

## To run on a test file

```shell
cd PhysicsTools/BParkingNANO/test/
cmsenv 
cmsRun test_BParkSequence_10215.py
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
