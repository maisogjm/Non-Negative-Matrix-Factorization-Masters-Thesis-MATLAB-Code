#!/bin/csh -f

# Obtain the date.
set theDate = `date`
set month   = $theDate[2]
set day     = $theDate[3]

setenv VERSION  NewBrunetOnly
setenv RUN      $argv[1]
setenv NUM_SIMS $argv[2]
setenv PROGDIR  /data/jmm97/PrototypeBigSimulation3/progress
setenv OUTDIR   /data/jmm97/PrototypeBigSimulation3/results

setenv OUTFILE   ${OUTDIR}/EstK.$VERSION.$HOST.$RUN.$NUM_SIMS.$month.$day.mat
setenv ERRFILE   ${PROGDIR}/Sim.${VERSION}.$HOST.$RUN.$NUM_SIMS.$month.$day.stderr.txt
setenv SEEDFILE  ${OUTDIR}/Seed.$VERSION.$HOST.$RUN.mat

if ( ! -e $ERRFILE ) then
    touch $ERRFILE
endif

/data/jmm97/PrototypeBigSimulation3/BrunetOnlyMultipleRuns >>& $ERRFILE

