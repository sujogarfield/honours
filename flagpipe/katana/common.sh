# Sourced by every flagpipe PBS job. Override any of these with qsub -v.
# Defaults live in your own scratch; point them at a shared folder to reuse
# someone else's env / eggNOG database.
FLAGPIPE_ENV=${FLAGPIPE_ENV:-/srv/scratch/$USER/flagpipe_env}
if [ -z "${CONDA_BASE:-}" ]; then
  if command -v conda >/dev/null 2>&1; then CONDA_BASE=$(conda info --base)
  else CONDA_BASE=/srv/scratch/$USER/miniforge3; fi
fi
export PATH="$FLAGPIPE_ENV/bin:$PATH"
FP="python flagpipe/flagpipe.py"
