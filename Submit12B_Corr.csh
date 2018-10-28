#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( SimWithCorrelatedDataNoBrunet1MNOP.sub SimWithCorrelatedDataNoBrunet1QRST.sub SimWithCorrelatedDataNoBrunet1UVWX.sub )
    condor_submit $s
end
