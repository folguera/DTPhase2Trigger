# DTPhase2Triger
```
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src/
cmsenv

git clone https://github.com/folguera/DTPhase2Trigger L1Trigger/DTPhase2Trigger

scramv1 b -j 5
cd L1Trigger/DTPhase2Trigger

python makeHoughTransform.py --input="path/to/your/DTNtuples.root" 
```
