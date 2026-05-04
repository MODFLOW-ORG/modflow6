cd "$GITHUB_WORKSPACE/petsc"
PETSC_PREFIX=$(cygpath -u "$GITHUB_WORKSPACE/petsc")/install
./configure \
    --prefix="$PETSC_PREFIX" \
    --with-debugging=0 \
    --with-shared-libraries=1 \
    --with-cc='cl' \
    --with-fc='ifort' \
    --with-cxx='cl' \
    --with-blaslapack-lib='-L/cygdrive/c/PROGRA~2/Intel/oneAPI/mkl/latest/lib mkl_intel_lp64_dll.lib mkl_sequential_dll.lib mkl_core_dll.lib' \
    --with-mpi-include='/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/include' \
    --with-mpi-lib='/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/lib/impi.lib' \
    --with-mpiexec='/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/bin/mpiexec -localonly' \
    --with-python-exec='/usr/bin/PYTHON~1.EXE' \
    --with-fortran-bindings='1'