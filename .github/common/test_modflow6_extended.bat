cd "%GITHUB_WORKSPACE%\modflow6\autotest"
pixi run autotest -m "%MARKERS%" -k "%FILTERS%" --netcdf --parallel