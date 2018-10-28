#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( FinalSimulationsABCD.NoCCC.sub FinalSimulationsEFGH.NoCCC.sub FinalSimulationsIJKL.NoCCC.sub )
    condor_submit $s
end
