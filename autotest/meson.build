# test-drive is built from source as a subproject so its module files
# always match the compiler building MF6. unit tests need gfortran.
# with unit_tests=auto, other compilers skip test-drive without
# downloading it; with unit_tests=enabled, failure is an error.
unit_tests = get_option('unit_tests').require(
  fc_id == 'gcc',
  error_message : 'Fortran unit tests require gfortran',
)
test_drive = dependency('', required : false)
test_drive_proj = subproject(
  'test-drive',
  required : unit_tests,
  # not our code, so don't fail on its warnings
  default_options : ['default_library=static', 'werror=false'],
)
if test_drive_proj.found()
  test_drive = test_drive_proj.get_variable('testdrive_dep')
endif
if test_drive.found()
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

  # the BMI suite runs the minimal simulation in a process of its own,
  # because a simulation cannot be initialized twice in one process
  test(
    'Bmi',
    tester,
    args: 'Bmi',
    workdir: meson.project_source_root() / '.mf6minsim',
    is_parallel: false,
  )
endif