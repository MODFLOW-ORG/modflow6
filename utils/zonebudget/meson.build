zonebudget_sources = files(
    'src' / 'budgetdata.f90',
    'src' / 'zbud6.f90',
    'src' / 'zone.f90',
    'src' / 'zoneoutput.f90',
)


zbud6name = 'zbud6'
message('ZONEBUDGET 6 executable name: ' + zbud6name + '\n')

executable(zbud6name,
           zonebudget_sources,
           link_with: mf6core, 
           install: true,)
