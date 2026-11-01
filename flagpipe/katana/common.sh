# Sourced by every flagpipe PBS job. Override any of these with qsub -v.
FLAGPIPE_ENV=${FLAGPIPE_ENV:-/srv/scratch/baker/flagpipe_env}
CONDA_BASE=${CONDA_BASE:-/srv/scratch/baker/miniforge3}
export PATH="$FLAGPIPE_ENV/bin:$PATH"
FP="python flagpipe/flagpipe.py"
