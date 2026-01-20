set FC=ifort
cd "%GITHUB_WORKSPACE%\modflow6"
pixi run setup -Dextended=true builddir
pixi run build builddir
pixi run test builddir
pixi run setup-mf5to6 builddir
pixi run build-mf5to6 builddir
