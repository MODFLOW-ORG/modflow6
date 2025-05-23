test_drive = dependency('test-drive', required : false)
if test_drive.found() and not fc_id.contains('intel')
  tests = [
    'ArrayHandlers',
    'DevFeature',
    'GeomUtil',
    'HashTable',
    'InputOutput',
    'KeyValueList',
    'List',
    'ListIterator',
    'MathUtil',
    'MemoryContainerIterator',
    'MemoryStore',
    'Message',
    'PtrHashTable',
    'Sim',
    'SVD',
    'SwfUtils',
    'TimeSelect',
    'TimeStepSelect',
    'UzfEtUtil'
  ]

  test_srcs = files(
    'tester.f90',
    'tester_utils.f90',
  )
  foreach t : tests
    test_srcs += files('Test@0@.f90'.format(t.underscorify()))
  endforeach

  tester = executable(
    'tester',
    sources: test_srcs,
    link_with: mf6core,
    dependencies: test_drive,
  )

  test('Test source modules', tester)

  foreach t : tests
    test(t, tester, args: t)
  endforeach
endif