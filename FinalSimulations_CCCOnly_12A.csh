#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( FinalSimulationsABCD.CCCOnly.sub FinalSimulationsEFGH.CCCOnly.sub FinalSimulationsIJKL.CCCOnly.sub )
    condor_submit $s
end
