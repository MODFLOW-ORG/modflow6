cd "$GITHUB_WORKSPACE/petsc"
PETSC_PREFIX=$(cygpath -u "$GITHUB_WORKSPACE/petsc")
make all
make check

# We need to make sure this folder doesn't exist otherwise the install target does nothing.
rm -rf "$PETSC_PREFIX/install"
make install