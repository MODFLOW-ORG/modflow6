bmi_sources = files(
    'bmi.f90',
    'mf6bmi.f90',
    'mf6bmiError.f90',
    'mf6bmiGrid.f90',
    'mf6bmiUtil.f90',
    'mf6xmi.F90',
    'mf6api.f90'
)

message('MODFLOW 6 shared library name:', 'lib' + buildname)

library(buildname, bmi_sources, link_with: mf6core, name_prefix: 'lib', install: true)
