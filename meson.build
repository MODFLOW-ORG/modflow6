project(
  'MODFLOW 6',
  'fortran',
  version: '6.7.0.dev2',
  license: 'CC0',
  meson_version: '>= 1.7.1',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'buildtype=release',
    'debug=false',
    'fortran_std=f2018',
])

if get_option('buildtype') == 'release'
  profile = 'release'
else
  profile = 'develop'
endif
message('Profile:', profile)

# parse compiler options
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
message('Compiler ID:', fc_id)
compile_args = []
link_args = []

# Command line options for gfortran
if fc_id == 'gcc'
  # General options
  compile_args += [
                   '-fall-intrinsics',
                   '-pedantic',
                   '-cpp',
                   '-Wcharacter-truncation',
                   '-Wno-unused-dummy-argument', # This makes problems with OOP
                   '-Wno-intrinsic-shadow',      # We shadow intrinsics with methods, which should be fine
                   '-Wno-maybe-uninitialized',   # "Uninitialized" flags produce false positives with allocatables
                   '-Wno-uninitialized',
                   ]
  
  # Define OS with gfortran for OS specific code
  # These are identical to pre-defined macros available with ifort
  system = build_machine.system()
  if system == 'linux'
    compile_args += '-D__linux__'
  elif system == 'darwin'
    compile_args += '-D__APPLE__'
  elif system == 'windows'
    compile_args += '-D_WIN32'
  endif
endif

# Command line options for ifort
if fc_id == 'intel-cl'
  # windows
  compile_args += ['/fpe:0',              # Activate all floating point exceptions
                   '/heap-arrays:0',
                   '/traceback',
                   '/fpp',                # Activate preprocessing
                   '/Qdiag-disable:7416', # f2008 warning
                   '/Qdiag-disable:7025', # f2008 warning
                   '/Qdiag-disable:5268', # Line too long
                   '/Qdiag-disable:10448',# ifort deprecation warning
                  ]
  link_args += ['/ignore:4217', # access through ddlimport might be inefficient
                '/ignore:4286'  # same as 4217, but more general
               ]
elif fc_id == 'intel'
  # linux and macOS
  compile_args += ['-fpe0',              # Activate all floating point exceptions
                   '-no-heap-arrays',
                   '-traceback',
                   '-diag-disable:7416',  # f2008 warning
                   '-diag-disable:7025',  # f2008 warning
                   '-diag-disable:5268',  # Line too long
                  ]
  link_args += '-static-intel'
  
# Command line options for ifx
elif fc_id == 'intel-llvm-cl'
  # windows
  compile_args += ['/fpe:0',              # Activate all floating point exceptions
                   '/heap-arrays:0',
                   '/traceback',
                   '/fpp',                # Activate preprocessing
                   '/Qdiag-disable:7416', # f2008 warning
                   '/Qdiag-disable:7025', # f2008 warning
                   '/Qdiag-disable:5268', # Line too long
                  ]
  link_args += ['/ignore:4217', # access through ddlimport might be inefficient
                '/ignore:4286'  # same as 4217, but more general
               ]

endif

# parallel build options
is_extended_build = get_option('extended')
is_parallel_build = get_option('parallel') or is_extended_build
is_netcdf_build = get_option('netcdf') or is_extended_build
is_cray = get_option('cray')
is_mpich = get_option('mpich')
if is_cray and build_machine.system() != 'linux'
  error('cray only supported on linux systems')
endif
if is_cray and not is_parallel_build
  is_extended_build = true
  is_parallel_build = true
  is_mpich = true
endif
if is_mpich
  if is_cray
    mpifort_name = 'mpichf90' 
  else
    mpifort_name = 'mpichfort'
  endif
endif

if is_extended_build
  message('Extended build:', is_extended_build)
elif is_parallel_build
  message('Parallel build:', is_parallel_build)
elif is_netcdf_build
  message('NetCDF build:', is_netcdf_build)
endif

# windows options for petsc
petsc_dir_rel = '..' /'petsc'
petsc_dir_abs = meson.project_source_root() / '..' /'petsc'
petsc_arch = 'arch-mswin-c-opt'
petsc_compiled_rel = petsc_dir_rel / petsc_arch
petsc_compiled_abs = petsc_dir_abs / petsc_arch

# windows options for netcdf
netcdf_dir_rel = '..' /'netcdf'
netcdf_dir_abs = meson.project_source_root() / '..' / 'netcdf'
netcdfc_build = 'netCDF4.9.2-NC4-64'
netcdff_build = 'netcdf-fortran-4.6.1' / 'build'
netcdfc_compiled_rel = netcdf_dir_rel / netcdfc_build
netcdfc_compiled_abs = netcdf_dir_abs / netcdfc_build
netcdff_compiled_rel = netcdf_dir_rel / netcdff_build
netcdff_compiled_abs = netcdf_dir_abs / netcdff_build

# on windows only with intel
if build_machine.system() == 'windows' and is_parallel_build
  if fc_id != 'intel-cl'
    error('Parallel build on Windows only with intel compiler. Terminating...')
  endif
endif

# lists for parallel dependencies and compiler arguments
dependencies = [ ]
extra_cmp_args = [ ]
extended_incdir = [ ]

# load petsc, mpi, and netcdf dependencies/libraries
if is_parallel_build or is_extended_build
  # find petsc
  if build_machine.system() != 'windows'
    petsc = dependency('PETSc', required : true)
  else
    # directly look for petsc
    petsc = fc.find_library('libpetsc', dirs: petsc_compiled_abs / 'lib', required : true)
  endif
  extra_cmp_args += [ '-D__WITH_PETSC__' ]
  dependencies += petsc
  with_petsc = true

  # find mpi
  if is_mpich
    mpifort = dependency(mpifort_name, required : true)
    dependencies += mpifort
  else
    mpi = dependency('mpi', language : 'fortran', required : true)
    dependencies += mpi
  endif
  extra_cmp_args += [ '-D__WITH_MPI__']
  with_mpi = true
else
  with_petsc = false
  with_mpi = false
endif

# find netcdf
with_netcdf = false
if is_netcdf_build
  netcdff = dependency('netcdf', language : 'fortran', static: false, required : true)
  if netcdff.found()
    with_netcdf = true
    extra_cmp_args += [ '-D__WITH_NETCDF__' ]
    dependencies += [ netcdff ]
    if build_machine.system() == 'windows'
      netcdf_incdir = include_directories([
        netcdfc_compiled_rel / 'include',
        netcdff_compiled_rel / 'fortran'
      ])
      extended_incdir += [
        '..' / netcdfc_compiled_rel / 'include',
        '..' / netcdff_compiled_rel / 'fortran'
      ]
      netcdff_libdir = netcdff_compiled_abs / 'fortran'
      lnetcdff = fc.find_library('netcdff', dirs: netcdff_libdir, static: false, required : true)
      dependencies += [ lnetcdff ]
    endif
  endif
endif

# GCC profile options need to be netcdf aware due to HDF5 issue
if fc_id == 'gcc'
  if with_netcdf
    # HDF5 1.14.3 invalid fpe trap issue: https://github.com/HDFGroup/hdf5/issues/3831
    if profile == 'release'
      compile_args += ['-ffpe-summary=overflow', '-ffpe-trap=overflow,zero']
    elif profile == 'develop'
      compile_args += ['-fcheck=all', '-ffpe-trap=overflow,zero']
    endif
  else
    if profile == 'release'
      compile_args += ['-ffpe-summary=overflow', '-ffpe-trap=overflow,zero,invalid']
    elif profile == 'develop'
      compile_args += ['-fcheck=all', '-ffpe-trap=overflow,zero,invalid']
    endif
  endif
endif

compile_args += extra_cmp_args

add_project_arguments(fc.get_supported_arguments(compile_args), language: 'fortran')
add_project_link_arguments(fc.get_supported_arguments(link_args), language: 'fortran')

if is_parallel_build and build_machine.system() == 'windows'
  message('Compiling PETSc Fortran modules')
  petsc_incdir = include_directories([petsc_dir_rel / 'include', petsc_compiled_rel / 'include'])
  extended_incdir += ['..' / petsc_dir_rel / 'include', '..' / petsc_compiled_rel / 'include']
  petsc_src = petsc_dir_abs / 'src'
  sources_petsc = [petsc_src / 'dm/f90-mod/petscdmdamod.F90',
                  petsc_src / 'dm/f90-mod/petscdmmod.F90',
                  petsc_src / 'dm/f90-mod/petscdmplexmod.F90',
                  petsc_src / 'dm/f90-mod/petscdmswarmmod.F90',
                  petsc_src / 'ksp/f90-mod/petsckspdefmod.F90',
                  petsc_src / 'ksp/f90-mod/petsckspmod.F90',
                  petsc_src / 'ksp/f90-mod/petscpcmod.F90',
                  petsc_src / 'mat/f90-mod/petscmatmod.F90',
                  petsc_src / 'snes/f90-mod/petscsnesmod.F90',
                  petsc_src / 'sys/f90-mod/petscsysmod.F90',
                  petsc_src / 'sys/mpiuni/f90-mod/mpiunimod.F90',
                  petsc_src / 'tao/f90-mod/petsctaomod.F90',
                  petsc_src / 'ts/f90-mod/petsctsmod.F90',
                  petsc_src / 'vec/f90-mod/petscvecmod.F90',]
  petsc_modules = static_library('petsc_modules', 
                                sources_petsc,
                                dependencies: dependencies,
                                include_directories: petsc_incdir)  
endif

# build mf6 and libmf6
buildname = get_option('buildname')
subdir('src')
subdir('srcbmi')

# build zbud6 and mf5to6 utility programs
subdir('utils/zonebudget')

# add autotest directory
fs = import('fs')
if fs.is_dir('autotest')
  subdir('autotest')
endif

# meson tests to evaluate installation success
testdir = meson.project_source_root() / '.mf6minsim'
if with_mpi
  mpiexec = find_program('mpiexec', required : false)
  if mpiexec.found()
    test('Parallel version command line test', mpiexec, args : ['-n', '2', mf6exe, '-v', '-p'], is_parallel : false)
    test('Parallel compiler command line test', mpiexec, args : ['-n', '2', mf6exe, '-c', '-p'], is_parallel : false)
    test('Serial simulation test', mf6exe, workdir : testdir, is_parallel : false)
    test('Parallel simulation test - 1 core', mpiexec, workdir : testdir, args : ['-n', '1', mf6exe, '-p'], is_parallel : false)
    test('Parallel simulation test - 2 cores', mpiexec, workdir : testdir, args : ['-n', '2', mf6exe, '-p'], is_parallel : false)
  endif
else
  test('Version command line test', mf6exe, args : ['-v',])
  test('Compiler command line test', mf6exe, args : ['-c',])
  test('Test installation help', mf6exe, args : ['-h',])
  test('Serial simulation test', mf6exe, workdir : testdir)
endif
