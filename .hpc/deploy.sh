#!/bin/bash

set -euxo pipefail

# fetch and checkout latest
ssh -l "$SSH_USERNAME" "$SSH_HOSTNAME" "cd $MF6_PROJ_ROOT && git fetch $GIT_REMOTE && git checkout $GIT_REMOTE/$GIT_BRANCH"
echo "Updated repository $MF6_PROJ_ROOT"
# submit a job to build mf6
jobid=$(ssh -l "$SSH_USERNAME" "$SSH_HOSTNAME" "sbatch --export=ALL,MESON_VERSION=${MESON_VERSION:-},PETSC_VERSION=${PETSC_VERSION:-} --account=$SLURM_ACCOUNT --reservation=$SLURM_RESERVATION --parsable -D $MF6_PROJ_ROOT $MF6_PROJ_ROOT/$BUILD_SCRIPT" | tail -n 1)
echo "Submitted build job $jobid"
# submit a job to update the mf6 module
jobid=$(ssh -l "$SSH_USERNAME" "$SSH_HOSTNAME" "sbatch --export=ALL,MF6_PREV_VERSION=$MF6_PREV_VERSION,MF6_PROJ_ROOT=$MF6_PROJ_ROOT,MODULES_PATH=$MODULES_PATH --account=$SLURM_ACCOUNT --reservation=$SLURM_RESERVATION --parsable -D $MF6_PROJ_ROOT -d afterok:$jobid $MF6_PROJ_ROOT/$MODULE_SCRIPT" | tail -n 1)
echo "Submitted module update job $jobid"
