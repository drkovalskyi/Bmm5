# Usage:
# sh job_starter.sh [Processor Type] [Job]

# Get temporary log
tmp_log=`mktemp --tmpdir tmpLogPPNA.XXXXXXX`

processor_name=$1
job_name=$2

final_log=${job_name/%.job/.log}

python job_starter.py $processor_name $job_name &> $tmp_log

mv $tmp_log $final_log
