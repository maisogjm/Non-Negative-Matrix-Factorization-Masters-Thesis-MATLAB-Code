#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( SimWithCorrelatedDataNoBrunet1ABCD.sub SimWithCorrelatedDataNoBrunet1EFGH.sub SimWithCorrelatedDataNoBrunet1IJKL.sub )
    condor_submit $s
end
