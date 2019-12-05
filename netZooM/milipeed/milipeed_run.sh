#!/bin/bash
# Working script for running milipeed in different scenarios.
# Edit 'milipeed_config.m' to set program parameters.

# Console running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('milipeed_config.m'); run('milipeed_run.m');"

# Foreground running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('milipeed_config.m'); run('milipeed_run.m'); quit;"

# Background running (./milipeed_run.sh &)
nohup matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('milipeed_config.m'); run('milipeed_run.m'); quit;" >& milipeed.`hostname`.log

# Email notification when done
echo "milipeed run on `hostname` has just finished: `date`." | mail -s "Task finished on `hostname`" -a milipeed.`hostname`.log `whoami`

# Testing
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('milipeed_config_test.m'); run('milipeed_run.m'); quit;"
#diff /tmp/milipeed.test.txt test_data/milipeed.txt

# Benchmark (5 repeats)
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('milipeed_config_test.m'); run('milipeed_run.m'); run('milipeed_run.m'); run('milipeed_run.m'); run('milipeed_run.m'); run('milipeed_run.m'); quit;" >& milipeed.test.`hostname`.log
