#!/bin/csh -f

# Check command line arguments.
if ( $#argv != 3 ) then
    echo "Need number of simulations, version, and run."
    exit 1
endif

# Set up environment variables to pass information to the MATLAB session.
setenv THE_NODE $HOSTNAME
setenv NUM_SIMS $argv[1] ; # Maybe this should be hard-coded within the MATLAB code itself.
setenv VERSION  $argv[2] ; # Maybe this should be hard-coded within the MATLAB code itself.
setenv RUN      $argv[3] ; # Maybe this should be hard-coded within the MATLAB code itself.
setenv OUTDIR   /home/jmm97/estimateKsims

# Obtain the date.
set theDate = `date`
set month   = $theDate[2]
set day     = $theDate[3]

# Determine output file name; pass this to the MATLAB session as an environment variable.
# Maybe this file name should be parsed together within the MATLAB code itself.
# The output file should contain not only simulation results, but also the
# starting random number seed, just in case we want to reproduce the results.
setenv OUTFILE /data/jmm97/PrototypeBigSimulation3/progress/EstK.$VERSION.$THE_NODE.$RUN.$NUM_SIMS.$month.$day.mat

# This is the name of the output file containing the randmo number seed for the next run.
# Early on, the MATLAB job will load this file into memory and set the random seed
# to the value contained therein, before generating random numbers.
setenv SEEDFILE /tmp/Seed.$THE_NODE.mat

# Invoke MATLAB simulation and processing.
#rajml -nojvm -nosplash -r PrototypeBigSimulation2
condor_submit PrototypeBigSimulation3.sub

# If the month is 'Jun' or if it is the first half of July, recursively submit a new job for tomorrow.
if ( 0 ) then
if ( ( $month == "Jun" ) || ( ( $month == "Jul" ) && ( $day <= 15 ) ) ) then

    # Schedule the next run.  Examples of using at:
    #     at 0815am Jan 24
    #     at 8:15am Jan 24
    #     at 3 pm Friday
    #     at now + 1 minute
    #     at now + 1 hour
    #     at now + 1 day
    at 3 am tomorrow << EOF
condor_submit ... *.sub
EOF
endif
endif

