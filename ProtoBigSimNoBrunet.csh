#!/bin/csh -f

# Obtain the date.
set theDate = `date`
set month   = $theDate[2]
set day     = $theDate[3]

setenv NUM_SIMS 20
setenv VERSION  NOBRUNET3
setenv RUN      $argv[1]
setenv PROGDIR  /data/jmm97/PrototypeBigSimulation3/progress
setenv OUTDIR   /data/jmm97/PrototypeBigSimulation3/results
setenv CCC_HIVAL 0

setenv OUTFILE   ${OUTDIR}/EstK.$VERSION.$HOST.$RUN.$NUM_SIMS.$month.$day.mat
setenv ERRFILE   ${PROGDIR}/Sim.$VERSION.$HOST.$RUN.$NUM_SIMS.$month.$day.stderr.txt
setenv SEEDFILE  ${OUTDIR}/Seed.$VERSION.$HOST.$RUN.mat

if ( ! -e $ERRFILE ) then
    touch $ERRFILE
endif

/data/jmm97/PrototypeBigSimulation3/ProtoBigSimNoBrunet >>& $ERRFILE

