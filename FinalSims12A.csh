#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( FinalSimulationsABCD.sub FinalSimulationsEFGH.sub FinalSimulationsIJKL.sub )
    condor_submit $s
end
