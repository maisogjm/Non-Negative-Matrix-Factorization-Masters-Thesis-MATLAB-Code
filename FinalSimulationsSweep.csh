#!/bin/csh -f

# Obtain the date.
set theDate = `date`
set month   = $theDate[2]
set day     = $theDate[3]

setenv VERSION  FinalSimulationsV7

setenv RUN       $argv[1]
setenv NUM_K     $argv[1]
setenv PROGDIR   /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations
setenv RESDIR    /data/jmm97/PrototypeBigSimulation3/results
setenv OUTDIR    ${RESDIR}/FinalSimulations
setenv SEEDDIR   ${RESDIR}/Seeds
#setenv OUTDIR   /data/jmm97/PrototypeBigSimulation3/results/FinalSimulations

setenv OUTFILE   ${OUTDIR}/EstK.$VERSION.$HOST.RUN.$NUM_K.$month.$day.mat
setenv ERRFILE   ${PROGDIR}/Sim.${VERSION}.$HOST.RUN.$NUM_K.$month.$day.stderr.txt
setenv SEEDFILE	 ${SEEDDIR}/Seed.$VERSION.$HOST.mat

if ( ! -e $ERRFILE ) then
    touch $ERRFILE
endif

/data/jmm97/PrototypeBigSimulation3/FinalSimulations >>& $ERRFILE

