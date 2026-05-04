cd "$GITHUB_WORKSPACE/petsc"
PETSC_PREFIX=$(cygpath -u "$GITHUB_WORKSPACE/petsc")
make all
make check
make PETSC_DIR=$PETSC_PREFIX PETSC_ARCH=arch-mswin-c-opt install