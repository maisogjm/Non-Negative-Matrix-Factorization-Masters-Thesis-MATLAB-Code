#!/bin/csh -f
cd /data/jmm97/PrototypeBigSimulation3
foreach s ( ProtoBigSimNoBrunetMNOP.sub ProtoBigSimNoBrunetQRST.sub ProtoBigSimNoBrunetUVWX.sub )
    condor_submit $s
end
