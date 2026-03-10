#!/bin/bash

LOGFILE="snakemake_$(date +%Y%m%d_%H%M%S).log"

nohup bash -c "
export PYTHONNOUSERSITE=True
snakemake --executor slurm \
          --default-resources slurm_account=carnegie_poc mem_mb=4000 runtime=240 \
          --rerun-incomplete \
	  --jobs 20
" > "$LOGFILE" 2>&1 &

# Save the process ID
echo $! > snakemake.pid

echo "Snakemake started with PID: $(cat snakemake.pid)"
echo "You can monitor progress with: tail -f $LOGFILE"
echo "You can check if it's still running with: ps -p \$(cat snakemake.pid)"
