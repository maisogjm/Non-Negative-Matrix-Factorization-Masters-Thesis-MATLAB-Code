#!/bin/csh -f

# Obtain the date.
set theDate = `date`
set month   = $theDate[2]
set day     = $theDate[3]

setenv VERSION  FinalSimulations_CCCOnlyV5

setenv RUN      $argv[1]
setenv PROGDIR  /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations
setenv OUTDIR   /data/jmm97/PrototypeBigSimulation3/results/FinalSimulations
setenv NUM_K    3

setenv OUTFILE   ${OUTDIR}/EstK.$VERSION.$HOST.$RUN.$NUM_K.$month.$day.mat
setenv ERRFILE   ${PROGDIR}/Sim.${VERSION}.$HOST.$RUN.$NUM_K.$month.$day.stderr.txt
setenv SEEDFILE  ${OUTDIR}/Seed.$VERSION.$HOST.mat

if ( ! -e $ERRFILE ) then
    touch $ERRFILE
endif

/data/jmm97/PrototypeBigSimulation3/FinalSimulations_CCCOnly >>& $ERRFILE

