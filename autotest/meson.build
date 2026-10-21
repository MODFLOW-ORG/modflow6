test_drive = dependency('test-drive', required : false)
if test_drive.found() and not fc_id.contains('intel')
  tests = [
    'ArrayHandlers',
    'BinaryFileReader',
    'Budget',
    'FeatureFlags',
    'GeomUtil',
    'GridFileReader',
    'HashTable',
    'InputOutput',
    'KeyValueList',
    'List',
    'ListIterator',
    'MathUtil',
    'MawNonVertical',
    'MawNur',
    'MemoryAlias',
    'MemoryContainerIterator',
    'MemoryStore',
    'Message',
    'PseudoInverse',
    'PtrHashTable',
    'Sim',
    'Smoothing',
    'SVD',
    'SwfUtils',
    'TimeSelect',
    'TimeStepSelect',
    'UzfCellGroup',
    'UzfEtUtil'
  ]

  test_srcs = files(
    'tester.f90',
    'tester_utils.f90',
    'TestBmi.f90',
    'TestXmi.f90',
  )
  test_srcs += bmi_sources
  foreach t : tests
    test_srcs += files('Test@0@.f90'.format(t.underscorify()))
  endforeach

  tester = executable(
    'tester',
    sources: test_srcs,
    link_with: mf6core,
    dependencies: [test_drive] + dependencies + extra_deps,
  )

  test('Test source modules', tester)

  foreach t : tests
    test(t, tester, args: t)
  endforeach

  # the BMI and XMI suites run the minimal simulation, each in a process of
  # its own, because a simulation cannot be initialized twice in one process
  foreach t : ['Bmi', 'Xmi']
    test(
      t,
      tester,
      args: t,
      workdir: meson.project_source_root() / '.mf6minsim',
      is_parallel: false,
    )
  endforeach
endif