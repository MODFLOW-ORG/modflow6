project(
  'mf5to6',
  'fortran',
  version: '6.7.0.dev1',
  license: 'CC0',
  meson_version: '>= 1.7.1',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'buildtype=release',
    'debug=false',
    'fortran_std=legacy',
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
                  ]
  link_args += ['/ignore:4217', # access through ddlimport might be inefficient
                '/ignore:4286'  # same as 4217, but more general
               ]
elif fc_id == 'intel'
  # linux and macOS
  compile_args += ['-fpe0',              # Activate all floating point exceptions
                   '-no-heap-arrays',
                   '-traceback',
                  ]
  link_args += '-static-intel'
  
# Command line options for ifx
elif fc_id == 'intel-llvm-cl'
  # windows
  compile_args += ['/fpe:0',              # Activate all floating point exceptions
                   '/heap-arrays:0',
                   '/traceback',
                   '/fpp',                # Activate preprocessing
                  ]
  link_args += ['/ignore:4217', # access through ddlimport might be inefficient
                '/ignore:4286'  # same as 4217, but more general
               ]

endif

# Suppress all warnings on mf5to6
if fc_id == 'gcc' or fc_id == 'intel'
  compile_args += '-w'
elif fc_id == 'intel-cl'
  compile_args += '/w'
endif 

add_project_arguments(fc.get_supported_arguments(compile_args), language: 'fortran')
add_project_link_arguments(fc.get_supported_arguments(link_args), language: 'fortran')

# build mf5to6 utility programs
subdir('src')
