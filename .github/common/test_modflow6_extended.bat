:: Copy Unidata netCDF/HDF5 DLLs into mf6 bin dir immediately before testing, which ensures the
:: MSVC-built Unidata DLLs are the last thing written to bin/ so they are linked at runtime.
copy /Y "%GITHUB_WORKSPACE%\netcdf\netCDF4.9.3-NC4-64\bin\hdf5.dll"    "%GITHUB_WORKSPACE%\modflow6\bin\"
copy /Y "%GITHUB_WORKSPACE%\netcdf\netCDF4.9.3-NC4-64\bin\hdf5_hl.dll" "%GITHUB_WORKSPACE%\modflow6\bin\"
copy /Y "%GITHUB_WORKSPACE%\netcdf\netCDF4.9.3-NC4-64\bin\netcdf.dll"  "%GITHUB_WORKSPACE%\modflow6\bin\"
copy /Y "%GITHUB_WORKSPACE%\netcdf\netCDF4.9.3-NC4-64\bin\zlib1.dll"   "%GITHUB_WORKSPACE%\modflow6\bin\"
copy /Y "%GITHUB_WORKSPACE%\netcdf\netCDF4.9.3-NC4-64\bin\libcurl.dll" "%GITHUB_WORKSPACE%\modflow6\bin\"
copy /Y "%GITHUB_WORKSPACE%\netcdf\netcdf-fortran-4.6.2\build\fortran\netcdff.dll" "%GITHUB_WORKSPACE%\modflow6\bin\"
cd "%GITHUB_WORKSPACE%\modflow6\autotest"
pixi run autotest -m "%MARKERS%" -k "%FILTERS%" --netcdf --parallel
