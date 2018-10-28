#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( ProtoBigSimNoBrunetABCD.sub ProtoBigSimNoBrunetEFGH.sub ProtoBigSimNoBrunetIJKL.sub )
    condor_submit $s
end
