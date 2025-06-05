## This is a CMake file, part of Unidata's netCDF Fortran package.
# Copyright 2012-2014, see the COPYRIGHT file for more information.

##################################
# Set Project Properties
##################################

#Minimum required CMake Version
cmake_minimum_required(VERSION 3.12.0)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

#Project Name for netCDF-Fortran
PROJECT (NC4F
LANGUAGES C Fortran
HOMEPAGE_URL "https://www.unidata.ucar.edu/software/netcdf/docs-fortran/"
DESCRIPTION "NetCDF is a set of software libraries and machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data."
)
SET(NC4F_CTEST_PROJECT_NAME "netcdf-fortran")

set(PACKAGE "${NC4F_CTEST_PROJECT_NAME}" CACHE STRING "")

#Project Version
SET(NC4F_VERSION_MAJOR 4)
SET(NC4F_VERSION_MINOR 6)

SET(NC4F_VERSION_PATCH 2)
SET(NC4F_VERSION_NOTE "")
SET(NC4F_VERSION ${NC4F_VERSION_MAJOR}.${NC4F_VERSION_MINOR}.${NC4F_VERSION_PATCH}${NC4F_VERSION_NOTE})
SET(VERSION ${NC4F_VERSION})
SET(NC4F_LIB_VERSION 7.1.0)
SET(NC4F_SO_VERSION 7)
SET(PACKAGE_VERSION ${VERSION})

# Set known dependencies.
SET(EXTRA_DEPS "")
set(TMP_BUILDNAME        "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_VERSION}-${CMAKE_SYSTEM_PROCESSOR}" CACHE STRING "Build name variable for CDash")

# For CMAKE_INSTALL_LIBDIR
INCLUDE(GNUInstallDirs)

###
# Allow for some customization of the buildname.
# This will make it easier to identify different builds,
# based on values passed from command line/shell scripts.
#
# For ctest scripts, we can use CTEST_BUILD_NAME.
###

SET(BUILDNAME_PREFIX "" CACHE STRING "")
SET(BUILDNAME_SUFFIX "" CACHE STRING "")

IF(BUILDNAME_PREFIX)
  SET(TMP_BUILDNAME "${BUILDNAME_PREFIX}-${TMP_BUILDNAME}")
ENDIF()

IF(BUILDNAME_SUFFIX)
  SET(TMP_BUILDNAME "${TMP_BUILDNAME}-${BUILDNAME_SUFFIX}")
ENDIF()

IF(NOT BUILDNAME)
  SET(BUILDNAME "${TMP_BUILDNAME}" CACHE STRING "Build name variable for CDash")
ENDIF()
###
# End BUILDNAME customization.
###

set(abs_top_builddir "${CMAKE_CURRENT_BINARY_DIR}")
set(abs_top_srcdir "${CMAKE_CURRENT_SOURCE_DIR}")

# Determine if DOXYGEN and DOT are available. These will be used
# when building the documentation.

FIND_PACKAGE(Doxygen)
FIND_PROGRAM(NC_DOT NAMES dot)

# A basic script used to convert m4 files
FIND_PROGRAM(NF_M4 NAMES m4 m4.exe)
IF(NF_M4)
  MESSAGE(STATUS "Found m4: ${NF_M4}")
  SET(HAVE_M4 TRUE)
ELSE()
  MESSAGE(STATUS "m4 not found.")
ENDIF()


# For CMAKE_INSTALL_LIBDIR
INCLUDE(GNUInstallDirs)

# Include this so we can check values in netcdf_meta.h.
INCLUDE(CheckCSourceCompiles)

# By default install Fortran modules into include directory
SET(Fortran_INSTALL_MODDIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE STRING "Fortran module install directory")

IF(MSVC)
  SET(GLOBAL PROPERTY USE_FOLDERS ON)
ENDIF()

# auto-configure style checks, other CMake modules.
INCLUDE (CheckLibraryExists)
INCLUDE (CheckIncludeFile)
INCLUDE (CheckTypeSize)
INCLUDE (CheckFunctionExists)
INCLUDE (CheckCXXSourceCompiles)
INCLUDE (TestBigEndian)
INCLUDE (CheckSymbolExists)
INCLUDE (GetPrerequisites)
INCLUDE (CheckFortranSourceCompiles)
INCLUDE (CheckFortranCompilerFlag)
FIND_PACKAGE(PkgConfig QUIET)

# Enable 'dist and distcheck'.
# File adapted from http://ensc.de/cmake/FindMakeDist.cmake
# TODO: this can't currently be found, need
# FIND_PACKAGE(MakeDist)
# End 'enable dist and distcheck'

# Set the build type.
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING "Choose the type of build, options are: None, Debug, Release."
    FORCE)
ENDIF()

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

################################
# End Project Properties
################################


################################
# Set CTest Properties
################################

ENABLE_TESTING()
INCLUDE(CTest)

# Copy the CTest customization file into binary directory, as required.
# TODO: need to provide suitable CTestCustom.cmake


# Set Memory test program for non-MSVC based builds.
# Assume valgrind for now.
IF(NOT MSVC)
  SET (CTEST_MEMORYCHECK_COMMAND valgrind CACHE STRING "")
ENDIF()

# Set variable to define the build type.
INCLUDE(GenerateExportHeader)

################################
# End CTest Properties
################################


################################
# Compiler Configuration
################################


# Default to shared libs on.
OPTION(BUILD_SHARED_LIBS "Configure netCDF as a shared library." ON)
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  SET(LIB_TYPE SHARED)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

# Supress unused variable and parameter warnings, for the time being,
# on GCC.
#
# Also, set some other default compiler flags.
IF(CMAKE_COMPILER_IS_GNUCC OR APPLE)
  string(APPEND CMAKE_C_FLAGS " -g -Wall -Wno-unused-variable -Wno-unused-parameter")
ENDIF()

IF(NOT ENABLE_COVERAGE_TESTS)
  string(APPEND CMAKE_C_FLAGS " -O2")
ENDIF()

IF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  IF(MSVC)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays:1")
  ENDIF()
ENDIF()

# Supress CRT Warnings.
# Only necessary for Windows
IF(MSVC)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
ENDIF()


################################
# End Compiler Configuration
################################


################################
# Define Utility Macros
################################

# Macro to append files to the EXTRA_DIST files.
SET(EXTRA_DIST "")
MACRO(ADD_EXTRA_DIST files)
  FOREACH(F IN LISTS files)
    list(APPEND EXTRA_DIST ${CMAKE_CURRENT_SOURCE_DIR}/${F})
  ENDFOREACH()
  SET(EXTRA_DIST ${EXTRA_DIST} PARENT_SCOPE)
ENDMACRO()

# A basic script used to convert m4 files
FIND_PROGRAM(NC4F_M4 NAMES m4)
MACRO(GEN_m4 filename)
  IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c
      COMMAND ${NC4F_M4}
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c
      VERBATIM
      )
  ENDIF()
ENDMACRO(GEN_m4)

# [LEFT OFF HERE] replacing macros with ones from netcdf-c macros using NC4F_prefix

# Binary tests, but ones which depend on value of 'TEMP_LARGE' being defined.
MACRO(add_bin_env_temp_large_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TEMP_LARGE=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()


# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Binary tests which are used by a script looking for a specific name.
MACRO(add_bin_test_no_prefix F)
  build_bin_test(${F})
  ADD_TEST(${F} ${EXECUTABLE_OUTPUT_PATH}/${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${F} PROPERTY FOLDER "tests/")
  ENDIF()
ENDMACRO()

MACRO(add_bin_test prefix F)

  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${F}.F)
    set(ext ".F")
  ELSE()
    set(ext ".f90")
  ENDIF()

  ADD_EXECUTABLE(${prefix}_${F} ${F}${ext})
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
  ADD_TEST(${prefix}_${F} ${EXECUTABLE_OUTPUT_PATH}/${prefix}_${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${prefix}_${F} PROPERTY FOLDER "tests/")
  ENDIF()
ENDMACRO()

# A script to print out information at the end of the configuration step.
MACRO(print_conf_summary)
  MESSAGE("")
  MESSAGE("")
  MESSAGE(STATUS "CMake Summary:")
  MESSAGE("")
  MESSAGE(STATUS "\tBuild Type:\t\t\t" ${CMAKE_BUILD_TYPE})
  MESSAGE(STATUS "\tBuilding Shared Libraries:\t" ${BUILD_SHARED_LIBS})
  MESSAGE(STATUS "\tBuilding netCDF-4:\t\t" ${ENABLE_NETCDF4})
  MESSAGE(STATUS "\tBuilding DAP Support:\t\t" ${ENABLE_DAP})
  MESSAGE(STATUS "\tLogging Support:\t\t" ${STATUS_LOGGING})
  MESSAGE(STATUS "\tEnabling CDF5:\t\t\t" ${ENABLE_CDF5})
  MESSAGE(STATUS "\tRelaxed check NC_EINVALCOORDS:\t" ${NF_RELAX_COORD_BOUND})
  MESSAGE(STATUS "\tBuilding Parallel IO:\t\t" ${BUILD_PARALLEL})
  MESSAGE(STATUS "\tUsing NetCDF4 Parallel IO:\t" ${NF_HAS_PARALLEL4})
  MESSAGE(STATUS "\tUsing PnetCDF Parallel IO:\t" ${NF_HAS_PNETCDF})
  MESSAGE(STATUS "\tTesting Parallel IO:\t\t" ${TEST_PARALLEL})

  IF(CMAKE_PREFIX_PATH)
    MESSAGE(STATUS "\tCMake Prefix Path:\t\t" "${CMAKE_PREFIX_PATH}")
  ENDIF()

  SET(ALL_DEP_LIBS "${NETCDF_C_LIBRARY};${EXTRA_DEPS}")
  MESSAGE(STATUS "\tLinking against:\t\t" "${ALL_DEP_LIBS}")
  FOREACH(_LIB IN LISTS EXTRA_DEPS)
    GET_FILENAME_COMPONENT(_LIB_NAME ${_LIB} NAME_WE)
    STRING(REGEX REPLACE "^lib" "" _NAME ${_LIB_NAME})
    LIST(APPEND NCF_LIBS "-l${_NAME}")
  ENDFOREACH()
  MESSAGE("")

ENDMACRO()
# Shell script Macro
MACRO(add_sh_test prefix F)
  IF(NOT MSVC)
    ADD_TEST(${prefix}_${F} bash "-c" "export srcdir=${CMAKE_CURRENT_SOURCE_DIR};export TOPSRCDIR=${CMAKE_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}/${F}.sh")
  ENDIF()

ENDMACRO()


# Macro for replacing '/MD' with '/MT'.
# Used only on Windows, /MD tells VS to use the shared
# CRT libs, MT tells VS to use the static CRT libs.
#
# Taken From:
#   http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F
#
MACRO(specify_static_crt_flag)

  SET(vars
    CMAKE_C_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(flag_var IN LISTS vars)
    IF(${flag_var} MATCHES "/MD")
      STRING(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    ENDIF()
  ENDFOREACH()

  FOREACH(flag_var IN LISTS vars)
    MESSAGE(STATUS " '${flag_var}': ${${flag_var}}")
  ENDFOREACH()
  MESSAGE(STATUS "")

ENDMACRO()

# A function used to create autotools-style 'yes/no' definitions.
# If a variable is set, it 'yes' is returned. Otherwise, 'no' is
# returned.
#
# Also creates a version of the ret_val prepended with 'NC',
# when feature is true, which is used to generate netcdf_meta.h.
FUNCTION(is_enabled feature ret_val)
  IF(${feature})
    SET(${ret_val} "yes" PARENT_SCOPE)
    SET("NC_${ret_val}" 1 PARENT_SCOPE)
  ELSE()
    SET(${ret_val} "no" PARENT_SCOPE)
    SET("NC_${ret_val}" 0 PARENT_SCOPE)
  ENDIF(${feature})
ENDFUNCTION()

# A function used to create autotools-style 'yes/no' definitions.
# If a variable is set, it 'yes' is returned. Otherwise, 'no' is
# returned.
#
# Also creates a version of the ret_val prepended with 'NC',
# when feature is true, which is used to generate netcdf_meta.h.
FUNCTION(is_disabled feature ret_val)
  IF(${feature})
    SET(${ret_val} "no" PARENT_SCOPE)
  ELSE()
    SET(${ret_val} "yes" PARENT_SCOPE)
    SET("NC_${ret_val}" 1 PARENT_SCOPE)
  ENDIF(${feature})
ENDFUNCTION()


################################
# End Macro Definitions
################################


#Add custom CMake Module
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/"
  CACHE INTERNAL "Location of our custom CMake modules.")

#####
# Some utility macros/scripts.

# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Build a binary used by a script, but don't make a test out of it.
MACRO(build_bin_test F ext)
  ADD_EXECUTABLE(${F} ${F}${ext})
  TARGET_LINK_LIBRARIES(${F} netcdff ${NETCDF_C_LIBRARY})
  IF(MSVC)
    SET_TARGET_PROPERTIES(${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
ENDMACRO()

OPTION(ENABLE_LARGE_FILE_TESTS "Run large file tests, which are slow and take lots of disk." OFF)
SET(LARGE_FILE_TESTS ${ENABLE_LARGE_FILE_TESTS})
OPTION(BUILD_BENCHMARKS "Run F90 I/O Benchmarks" OFF)
OPTION(TEST_WITH_VALGRIND "Run extra tests with valgrind" OFF)
OPTION(ENABLE_PARALLEL_TESTS "Run parallel I/O tests for F90 and F77" OFF)

SET(TEST_PARALLEL OFF)
IF(ENABLE_PARALLEL_TESTS)
  SET(TEST_PARALLEL ON)
ENDIF()

###
# Seek out netcdf library file and include path.
###

IF(NOT netCDF_LIBRARIES AND NOT netCDF_INCLUDE_DIR)
  FIND_PACKAGE(netCDF QUIET)
ELSE()
  SET(netCDF_FOUND TRUE)
ENDIF()

IF (netCDF_FOUND)
  SET(NETCDF_C_LIBRARY ${netCDF_LIBRARIES})
  SET(NETCDF_C_INCLUDE_DIR ${netCDF_INCLUDE_DIR})
  include_directories(${NETCDF_C_INCLUDE_DIR})
  message(STATUS "Found netCDF CMake package: ${NETCDF_C_LIBRARY}")
ELSE()
  # netCDF not properly packaged. Try to find it manually.
  FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf libnetcdf)
  IF(NOT NETCDF_C_LIBRARY)
    MESSAGE(FATAL_ERROR "libnetcdf not found. Please reinstall and try again.")
  ELSE()
    FIND_PATH(NETCDF_C_INCLUDE_DIR "netcdf.h")
    MESSAGE(STATUS "Found netcdf: ${NETCDF_C_LIBRARY}")
    IF(NOT NETCDF_C_INCLUDE_DIR)
      MESSAGE(FATAL_ERROR "Directory containing netcdf.h cannot be found. Please reinstall and try again.")
    ELSE()
      GET_FILENAME_COMPONENT(netCDF_LIB_DIR "${NETCDF_C_LIBRARY}" REALPATH)
    ENDIF()
  ENDIF()
ENDIF()

IF (NETCDF_C_INCLUDE_DIR)
   file(READ "${NETCDF_C_INCLUDE_DIR}/netcdf_meta.h" header)
   string(REGEX MATCH "#define NC_HAS_HDF5 *[01]" macrodef "${header}")
   string(REGEX MATCH "[01]" USE_NETCDF4 "${macrodef}")
   IF (USE_NETCDF4)
      MESSAGE(STATUS "Whether NetCDF-C built with HDF5 enabled: yes")
      #FIND_PACKAGE(HDF5 COMPONENTS C HL REQUIRED)
   ELSE()
      MESSAGE(STATUS "Whether NetCDF-C built with HDF5 enabled: no")
   ENDIF()
ENDIF()

if (NOT TARGET netCDF::netcdf)
  add_library(netCDF::netcdf UNKNOWN IMPORTED)
  set_target_properties(netCDF::netcdf PROPERTIES
    IMPORTED_LOCATION "${NETCDF_C_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_C_INCLUDE_DIR}"
  )
endif()

###
# End 'Seek out netcdf library file and include path.'
###

###
# Check for various header files
CHECK_INCLUDE_FILE("alloca.h"  HAVE_ALLOCA_H)
CHECK_INCLUDE_FILE("malloc.h"    HAVE_MALLOC_H)
CHECK_INCLUDE_FILE("ctype.h"   HAVE_CTYPE_H)
CHECK_INCLUDE_FILE("dirent.h"  HAVE_DIRENT_H)
CHECK_INCLUDE_FILE("dlfcn.h"   HAVE_DLFCN_H)
CHECK_INCLUDE_FILE("errno.h"   HAVE_ERRNO_H)
CHECK_INCLUDE_FILE("fcntl.h"   HAVE_FCNTL_H)
CHECK_INCLUDE_FILE("getopt.h"  HAVE_GETOPT_H)
CHECK_INCLUDE_FILE("stdbool.h"   HAVE_STDBOOL_H)
CHECK_INCLUDE_FILE("locale.h"  HAVE_LOCALE_H)
CHECK_INCLUDE_FILE("stdint.h"  HAVE_STDINT_H)
CHECK_INCLUDE_FILE("stdio.h"   HAVE_STDIO_H)
CHECK_INCLUDE_FILE("stdlib.h"  HAVE_STDLIB_H)
CHECK_INCLUDE_FILE("stdarg.h"    HAVE_STDARG_H)
CHECK_INCLUDE_FILE("strings.h"   HAVE_STRINGS_H)
CHECK_INCLUDE_FILE("signal.h"    HAVE_SIGNAL_H)
CHECK_INCLUDE_FILE("sys/dir.h"   HAVE_SYS_DIR_H)
CHECK_INCLUDE_FILE("sys/ndir.h"  HAVE_SYS_NDIR_H)
CHECK_INCLUDE_FILE("sys/param.h" HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE("sys/stat.h"  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE("sys/time.h"  HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE("sys/wait.h"  HAVE_SYS_WAIT_H)
CHECK_INCLUDE_FILE("sys/resource.h" HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE("fcntl.h"  HAVE_FCNTL_H)
CHECK_INCLUDE_FILE("inttypes.h"  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE("pstdint.h"  HAVE_PSTDINT_H)
CHECK_INCLUDE_FILE("endian.h" HAVE_ENDIAN_H)
CHECK_INCLUDE_FILE("BaseTsd.h"  HAVE_BASETSD_H)
CHECK_INCLUDE_FILE("stddef.h"   HAVE_STDDEF_H)
CHECK_INCLUDE_FILE("memory.h" HAVE_MEMORY_H)
SET(CMAKE_REQUIRED_INCLUDES ${NETCDF_C_INCLUDE_DIR})
CHECK_INCLUDE_FILE("netcdf.h" HAVE_NETCDF_H)


###
# End Check for various system header files.
###

add_compile_definitions(NCSHORT_T=integer*2)
add_compile_definitions(NCBYTE_T=byte)
add_compile_definitions(NF_RELAX_COORD_BOUND)
add_compile_definitions(NF_HAS_PARALLEL4)
add_compile_definitions(NF_HAS_PNETCDF)
#add_compile_definitions(ENABLE_CDF5)
add_compile_definitions(NF_INT1_IS_C_SIGNED_CHAR)
#add_compile_definitions(NF_INT1_IS_C_SHORT)
#add_compile_definitions(NF_INT1_IS_C_INT)
#add_compile_definitions(NF_INT1_IS_C_LONG)
add_compile_definitions(NF_INT2_IS_C_SHORT)
#add_compile_definitions(NF_INT2_IS_C_INT)
#add_compile_definitions(NF_INT2_IS_C_LONG)
add_compile_definitions(NF_INT_IS_C_INT)
#add_compile_definitions(NF_INT_IS_C_LONG)
add_compile_definitions(NF_INT8_IS_C_LONG_LONG)
add_compile_definitions(NF_REAL_IS_C_FLOAT)
#add_compile_definitions(NF_REAL_IS_C_DOUBLE)
add_compile_definitions(NF_DOUBLEPRECISION_IS_C_DOUBLE)
#add_compile_definitions(NF_DOUBLEPRECISION_IS_C_FLOAT)

###
# Type checks
###
CHECK_TYPE_SIZE("void*"     SIZEOF_VOIDSTAR)
CHECK_TYPE_SIZE("char"      SIZEOF_CHAR)
CHECK_TYPE_SIZE("double"    SIZEOF_DOUBLE)
CHECK_TYPE_SIZE("float"     SIZEOF_FLOAT)
CHECK_TYPE_SIZE("int"       SIZEOF_INT)
CHECK_TYPE_SIZE("long"      SIZEOF_LONG)
CHECK_TYPE_SIZE("long long" SIZEOF_LONG_LONG)
CHECK_TYPE_SIZE("off_t"     SIZEOF_OFF_T)
CHECK_TYPE_SIZE("off64_t"   SIZEOF_OFF64_T)
CHECK_TYPE_SIZE("short"     SIZEOF_SHORT)
CHECK_TYPE_SIZE("size_t"    SIZEOF_SIZE_T)
CHECK_TYPE_SIZE("ssize_t"   SIZEOF_SSIZE_T)
# __int64 is used on Windows for large file support.
CHECK_TYPE_SIZE("__int64"   SIZEOF___INT_64)
CHECK_TYPE_SIZE("uchar"     SIZEOF_UCHAR)
CHECK_TYPE_SIZE("int64_t"   SIZEOF_INT64_T)
CHECK_TYPE_SIZE("uint64_t"  SIZEOF_UINT64_T)

# On windows systems, we redefine off_t as __int64
# to enable LFS. This is true on 32 and 64 bit system.s
# We must redefine SIZEOF_OFF_T to match.
IF(MSVC AND SIZEOF___INT_64)
  SET(SIZEOF_OFF_T  ${SIZEOF___INT_64})
ENDIF()

IF(SIZEOF_SSIZE_T)
  SET(HAVE_SSIZE_T 1)
ELSE()
  CHECK_TYPE_SIZE("SSIZE_T" SIZEOF_SSIZE_T)
  IF(SIZEOF_SSIZE_T)
    SET(HAVE_SSIZE_T 1)
  ENDIF()
ENDIF()


###
# End type checks
###

###
# TODO: Replace these with automatically-determined
# checks, like in configure.ac.
###

##
# Check to see if compiler supports
# -fallow-argument-mismatch-flag introduced in gcc 10.
#
# Also check for NAG option -mismatch_all
#
# See https://github.com/Unidata/netcdf-fortran/issues/212
# See https://github.com/Unidata/netcdf-fortran/issues/218
##
check_Fortran_compiler_flag("-fallow-argument-mismatch" COMPILER_HAS_ALLOW_ARGUMENT_MISMATCH)

IF(COMPILER_HAS_ALLOW_ARGUMENT_MISMATCH)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
ENDIF(COMPILER_HAS_ALLOW_ARGUMENT_MISMATCH)

check_Fortran_compiler_flag("-mismatch_all" COMPILER_HAS_MISMATCH_ALL)

IF(COMPILER_HAS_MISMATCH_ALL)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mismatch_all")
ENDIF(COMPILER_HAS_MISMATCH_ALL)

CHECK_FUNCTION_EXISTS(alloca HAVE_ALLOCA)

CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nccreate "" USE_NETCDF_V2)
CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_set_log_level "" USE_LOGGING)
CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} oc_open "" BUILD_DAP)


###
# Check to see if szip write capability is present in netcdf-c.
###
SET(CMAKE_REQUIRED_INCLUDES ${NETCDF_C_INCLUDE_DIR})
CHECK_C_SOURCE_COMPILES("
#include <netcdf_meta.h>
#if !NC_HAS_SZIP_WRITE
      choke me
#endif
int main() {return 0;}" HAVE_SZIP_WRITE)

###
# Check to see if quantize capability is present in netcdf-c.
###
SET(CMAKE_REQUIRED_INCLUDES ${NETCDF_C_INCLUDE_DIR})
CHECK_C_SOURCE_COMPILES("
#include <netcdf_meta.h>
#if !NC_HAS_QUANTIZE
      choke me
#endif
int main() {return 0;}" HAVE_QUANTIZE)
if (HAVE_QUANTIZE)
   add_compile_definitions(NF_HAS_QUANTIZE)
endif()

CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_def_var_szip "" HAVE_DEF_VAR_SZIP)
if (NOT HAVE_DEF_VAR_SZIP)
  message(FATAL_ERROR "netcdf-c version 4.7.4 or greater is required."
                      "If netCDF-C 4.7.4 is present but is built statically, you *must* specify all dependencies"
                      "by setting the LDFLAGS environmental variable.\n"
                      "Example:\n"
                      "    $ export LDFLAGS=$(nc-config --libs)"
                      )
endif()

###
# ZStd Stanza
###

###
# Check to see if zstandard compression capability is present in netcdf-c.
###

OPTION(DISABLE_ZSTANDARD_PLUGIN "Disable ZStandard plugin detection and functionality, even if libnetcdf was compiled with plugin support" OFF)

IF(NOT DISABLE_ZSTANDARD_PLUGIN)
  SET(CMAKE_REQUIRED_INCLUDES ${NETCDF_C_INCLUDE_DIR})
  CHECK_C_SOURCE_COMPILES("
  #include <netcdf_meta.h>
  #if !NC_HAS_ZSTD
        choke me
  #endif
  int main() {return 0;}" NC_HAVE_ZSTD)

  IF(NC_HAVE_ZSTD)
    IF(NOT DEFINED ENV{HDF5_PLUGIN_PATH})
      MESSAGE(WARNING "libnetcdf was built with zstd support, but HDF5_PLUGIN_PATH is not set.  Please set the HDF5_PLUGIN_PATH environmental variable if zstandard tests fail, or pass -DDISABLE_ZSTANDARD_PLUGIN=TRUE when running cmake.")
      SET(NF_HAVE_ZSTD FALSE)
      SET(NF_HDF5_PLUGIN_PATH "")
    ELSE()
      MESSAGE(STATUS "HDF5_PLUGIN_PATH: $ENV{HDF5_PLUGIN_PATH}")
      SET(NF_HAVE_ZSTD TRUE)
      SET(NF_HDF5_PLUGIN_PATH "(HDF5_PLUGIN_PATH: $ENV{HDF5_PLUGIN_PATH}")
    ENDIF()
  ENDIF(NC_HAVE_ZSTD)
ELSE(NOT DISABLE_ZSTANDARD_PLUGIN)
  SET(NF_HAVE_ZSTD FALSE)
  SET(NF_HDF5_PLUGIN_PATH "")
ENDIF(NOT DISABLE_ZSTANDARD_PLUGIN)
###
# End ZStd Stanza
###


SET(HAS_NC4 FALSE)
IF(USE_NETCDF4)
  SET(NC_BUILD_V4 TRUE)
  SET(HAVE_NC_DEF_OPAQUE TRUE)
  SET(HAS_NC4 TRUE)
ENDIF(USE_NETCDF4)

SET(HAS_NC2 FALSE)
IF(USE_NETCDF_V2)
  SET(HAVE_NCCREATE TRUE)
  SET(USE_NC2 TRUE)
  SET(HAS_NC2 TRUE)
ENDIF(USE_NETCDF_V2)

SET(STATUS_LOGGING FALSE)
IF(USE_LOGGING)
  SET(HAVE_NC_SET_LOG_LEVEL TRUE)
  SET(STATUS_LOGGING ON)
ENDIF(USE_LOGGING)

SET(ENABLE_DAP OFF)
SET(HAS_DAP FALSE)
IF(BUILD_DAP)
  SET(HAVE_OC_OPEN TRUE)
  SET(ENABLE_DAP ON)
  SET(HAS_DAP TRUE)
ENDIF(BUILD_DAP)

##
# Set some status messages for the summary.
##


##
# End status messages/variables.
##

SET(ENABLE_NETCDF4 OFF)
IF(USE_NETCDF4)
  SET(ENABLE_NETCDF4 ON)
  add_compile_definitions(USE_NETCDF4)
  SET(ENABLE_NETCDF4 ON CACHE BOOL "")
ENDIF()

IF (UNIX AND ${CMAKE_SIZEOF_VOID_P} MATCHES "8")
  add_compile_definitions(LONGLONG_IS_LONG)
  # for easier debugging of cfortran.h
  IF (APPLE)
    IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmacro-backtrace-limit=0")
    ENDIF()
  ENDIF(APPLE)
ENDIF()

# Determine C library setting for NC_RELAX_COORD_BOUND
try_compile(NF_RELAX_COORD_BOUND ${CMAKE_CURRENT_BINARY_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/check_relax_coord_bound.c"
  LINK_LIBRARIES netCDF::netcdf
  )

# Determine C library setting for NC_HAS_CDF5
try_compile(ENABLE_CDF5 ${CMAKE_CURRENT_BINARY_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/check_cdf5.c"
  LINK_LIBRARIES netCDF::netcdf
  )

# Determine C library setting for NC_HAS_PNETCDF
try_compile(NF_HAS_PNETCDF ${CMAKE_CURRENT_BINARY_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/check_pnetcdf.c"
  LINK_LIBRARIES netCDF::netcdf
  )

# Determine C library setting for NC_HAS_PARALLEL4
try_compile(NF_HAS_PARALLEL4 ${CMAKE_CURRENT_BINARY_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/check_parallel4.c"
  LINK_LIBRARIES netCDF::netcdf
  )

IF(NF_HAS_PARALLEL4)
  add_compile_definitions(NF_HAS_PARALLEL4)
ENDIF()
IF(NF_HAS_PNETCDF)
  add_compile_definitions(NF_HAS_PNETCDF)
ENDIF()

IF(NF_HAS_PNETCDF OR NF_HAS_PARALLEL4)
  FIND_PACKAGE(MPI REQUIRED)
  IF(NOT MPI_Fortran_FOUND)
    MESSAGE(FATAL_ERROR
     " -----------------------------------------------------------------------\n"
     "  The NetCDF C library is built with parallel I/O feature enabled, but\n"
     "  the Fortran compiler supplied in this configure command\n"
     "  does not support MPI-IO. Please use one that does. If parallel I/O\n"
     "  feature is not desired, please use a NetCDF C library with parallel\n"
     "  I/O feature disabled. Abort.\n"
     " -----------------------------------------------------------------------")
  ENDIF(NOT MPI_Fortran_FOUND)
  SET(BUILD_PARALLEL ON)
  SET(HAVE_NC_USE_PARALLEL_ENABLED TRUE)
ELSE(NF_HAS_PNETCDF OR NF_HAS_PARALLEL4)
  SET(BUILD_PARALLEL OFF)
  SET(HAVE_NC_USE_PARALLEL_ENABLED FALSE)
  SET(TEST_PARALLEL OFF)
ENDIF(NF_HAS_PNETCDF OR NF_HAS_PARALLEL4)

# Determine C/Fortran pointer compatibility.
try_compile(COMPILE_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/test_c_ptrdiff_t.f90"
  )
IF (COMPILE_SUCCESS)
  add_compile_definitions(HAVE_TS29113_SUPPORT)
  SET(HAVE_TS29113_SUPPORT TRUE)
ELSE ()
  # If the compile fails, PTRDIFF_T is not defined.
  # Get sizeof of ptrdiff_t.
  SET(EXEC_NAME "${CMAKE_CURRENT_BINARY_DIR}/sizeof_ptrdiff_t")
  try_compile(COMPILE_SUCCESS2 ${CMAKE_CURRENT_BINARY_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/sizeof_ptrdiff_t.c"
    COPY_FILE ${EXEC_NAME}
    )
  IF (COMPILE_SUCCESS2)
    EXECUTE_PROCESS(
      COMMAND ${EXEC_NAME}
      OUTPUT_VARIABLE SIZEOF_PTRDIFF_T
      RESULT_VARIABLE SIZEOF_RESULT
      )
    IF (SIZEOF_RESULT)
      MESSAGE (FATAL_ERROR "UNABLE TO DETERMINE SIZEOF PTRDIFF_T")
    ELSE (SIZEOF_RESULT)
      add_compile_definitions(SIZEOF_PTRDIFF_T=${SIZEOF_PTRDIFF_T})
    ENDIF (SIZEOF_RESULT)
  ELSE()
    MESSAGE(FATAL_ERROR "Unable to compile ptrdiff")
  ENDIF()
endif()

IF(USE_NETCDF_V2)
  SET(BUILD_V2 ON CACHE BOOL "")
ENDIF()

OPTION(ENABLE_FORTRAN_TYPE_CHECKS
  "Determine Fortran types corresponding to netCDF types" ON)

# Toggle whether or not to run tests.
OPTION(ENABLE_TESTS "Enable netcdf-fortran tests." ON)

# Enable filter testing
OPTION(ENABLE_FILTER_TEST "Enable filter test." OFF)

# Set the default fortran builds; default is to build f03
SET(BUILD_F90 "ON")
SET(BUILD_V2 "ON")
SET(BUILD_F03 "ON")

IF(BUILD_F90)
  SET(HAS_F90 TRUE)
ELSE()
  SET(HAS_F90 FALSE)
ENDIF(BUILD_F90)

IF(BUILD_F03)
  SET(HAS_F03 TRUE)
ELSE()
  SET(HAS_F03 FALSE)
ENDIF(BUILD_F03)

IF(ENABLE_FORTRAN_TYPE_CHECKS)
  # Do tests to determine which Fortran types correspond to NCBYTE, NCSHORT, ...
  # The configure file got these by testing an F77 program, invoking
  # UD_FORTRAN_TYPES defined in acinclude.m4.

  INCLUDE(CMakeExtras/MatchNetCDFFortranTypes.cmake)
ELSE()
  # Set Fortran types to default.
  SET(NCBYTE_T "byte")
  SET(NCBYTE_T "byte")
  SET(NCSHORT_T "integer*2")
  SET(NF_INT1_T "byte")
  SET(NF_INT2_T "integer*2")
  SET(NF_INT8_T "integer*8")
  # Default is for following to be true.  TODO: test instead
  SET(NF_INT1_IS_C_SIGNED_CHAR "1")
  SET(NF_INT2_IS_C_SHORT "1")
  SET(NF_INT_IS_C_INT "1")
  SET(NF_INT8_IS_C_LONG_LONG "1")
  SET(NF_REAL_IS_C_FLOAT "1")
  SET(NF_DOUBLEPRECISION_IS_C_DOUBLE "1")
ENDIF(ENABLE_FORTRAN_TYPE_CHECKS)

# Create libsrc/nfconfig.inc (don't need nfconfig1.inc temporary file used with autoconf)
# configure_file("${NC4F_SOURCE_DIR}/libsrc/nfconfig.in.cmake"
#   "${NC4F_BINARY_DIR}/libsrc/nfconfig.inc")
# # Create config.h file
# configure_file("${NC4F_SOURCE_DIR}/config.h.cmake.in"
#   "${NC4F_BINARY_DIR}/config.h")

# make sure previous two files get cleaned up...
#SET_DIRECTORY_PROPERTIES (DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${netcdf-fortran_BINARY_DIR}/config.h ${netcdf-fortran_BINARY_DIR}/nf-config)

##
# Configuration for post-install RPath
# Adapted from http://www.cmake.org/Wiki/CMake_RPATH_handling
##
IF(NOT MSVC)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  IF(APPLE)
    SET(CMAKE_MACOSX_RPATH ON)
  ENDIF()

  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  ENDIF("${isSystemDir}" STREQUAL "-1")
ENDIF()

##
# End configuration for post-install RPath
##

###
# Test-related options
###

# Set CDash-related options.
#SET(NC4F_CTEST_DROP_SITE "my.cdash.org" CACHE STRING "Dashboard location for CTest-based testing purposes.")
SET(NC4F_CTEST_DROP_SITE "cdash.unidata.ucar.edu" CACHE STRING "Dashboard location for CTest-based testing purposes.")
SET(NC4F_CTEST_DROP_LOC_PREFIX "" CACHE STRING "Prefix for Dashboard location on remote server when using CTest-based testing.")

# Generate CTestConfig file from template.
CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/CTestConfig.cmake.in"
  "${CMAKE_SOURCE_DIR}/CTestConfig.cmake"
  @ONLY)

# Fortran compiler flag for the (module) include directory. Used for generating nf-config.
SET(MOD_FLAG -I)

###https://github.com/Unidata/netcdf-fortran/issues/321
# Allow the user to specify libraries
# to link against, similar to automakes 'LIBS' variable.
###
SET(NC_EXTRA_DEPS "" CACHE STRING "Additional libraries to link against.")
IF(NC_EXTRA_DEPS)
  STRING(REPLACE " " ";" DEPS_LIST ${NC_EXTRA_DEPS})
  FOREACH(_DEP IN LISTS DEPS_LIST)
    STRING(REGEX REPLACE "^-l" "" _LIB ${_DEP})
    FIND_LIBRARY("${_LIB}_DEP" NAMES "${_LIB}" "lib${_LIB}")
    # MESSAGE(${${_LIB}_DEP})
    IF(${_LIB}_DEP STREQUAL "")
      MESSAGE(FATAL_ERROR "Error finding ${_LIB}.")
    ELSE()
      MESSAGE(STATUS "Found ${_LIB}: ${${_LIB}_DEP}")
    ENDIF()
    SET(EXTRA_DEPS ${EXTRA_DEPS} "${${_LIB}_DEP}")
  ENDFOREACH()
  MESSAGE("Extra deps: ${EXTRA_DEPS}")
  LIST(REMOVE_DUPLICATES EXTRA_DEPS)
ENDIF()
###
# End user-specified dependent libraries.
###

###
# Option to enable examples.
###
OPTION(BUILD_EXAMPLES "Enable compilation of examples." ON)

# Determine whether or not to generate documentation.
OPTION(ENABLE_DOXYGEN "Enable generation of doxygen-based documentation." OFF)
IF(ENABLE_DOXYGEN)
  IF(DOXYGEN_FOUND)
    # Offer the option to build internal documentation.
    OPTION(ENABLE_INTERNAL_DOCS "Build internal documentation. This is of interest to developers only." OFF)
    IF(ENABLE_INTERNAL_DOCS)
      SET(BUILD_INTERNAL_DOCS YES CACHE STRING "")
    ELSE()
      SET(BUILD_INTERNAL_DOCS NO CACHE STRING "")
    ENDIF()

    # Option to turn on the TODO list in the doxygen-generated documentation.
    OPTION(ENABLE_DOXYGEN_TASKS "Turn on test, todo, bug lists in documentation. This is of interest to developers only." OFF)
    IF(ENABLE_DOXYGEN_TASKS)
      SET(SHOW_DOXYGEN_TAG_LIST YES CACHE STRING "")
    ELSE(ENABLE_DOXYGEN_TASKS)
      SET(SHOW_DOXYGEN_TODO_LIST NO CACHE STRING "")
    ENDIF(ENABLE_DOXYGEN_TASKS)

    # Specify whether or not 'dot' was found on the system path.
    IF(NC_DOT)
      SET(HAVE_DOT YES CACHE STRING "")
    ELSE(NC_DOT)
      SET(HAVE_DOT NO CACHE_STRING "")
    ENDIF(NC_DOT)

  ELSE(DOXYGEN_FOUND)
    MESSAGE(STATUS "Unable to build internal documentation.  Doxygen does not appear to be on your executable path. Install doxygen and configure the project again.")
    SET(ENABLE_DOXYGEN OFF)
  ENDIF(DOXYGEN_FOUND)
ENDIF()

###
# Set global linker flags
###
IF(CMAKE_LINK_FLAGS)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_LINK_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_LINK_FLAGS}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CMAKE_LINK_FLAGS}")
ENDIF()

# End autotools-style checks for config.h

# Recurse into other subdirectories.
ADD_SUBDIRECTORY(fortran)

ADD_SUBDIRECTORY(libsrc)

IF(ENABLE_TESTS)
  ADD_SUBDIRECTORY(nf_test)
  ADD_SUBDIRECTORY(nf03_test)
  ADD_SUBDIRECTORY(nf_test4)
  ADD_SUBDIRECTORY(nf03_test4)
ENDIF()

# Build the Fortran examples
IF(BUILD_EXAMPLES)
  ADD_SUBDIRECTORY(examples)
ENDIF(BUILD_EXAMPLES)
ADD_SUBDIRECTORY(docs)
print_conf_summary()

#####
# Configure and print the libnetcdff.settings file.
#####

# Set variables to mirror those used by autoconf.
# This way we don't need to maintain two separate template
# files.
SET(host_cpu "${CMAKE_SYSTEM_PROCESSOR}")
SET(host_vendor "${CMAKE_SYSTEM_NAME}")
SET(host_os "${CMAKE_SYSTEM_VERSION}")
SET(abs_top_builddir "${CMAKE_CURRENT_BINARY_DIR}")
SET(abs_top_srcdir "${CMAKE_CURRENT_SOURCE_DIR}")

SET(FC_VERSION "${CMAKE_Fortran_COMPILER}")

# Build *FLAGS for libnetcdf.settings.
SET(FFLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
SET(LDFLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")
SET(CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
SET(CPPFLAGS "${CMAKE_CPP_FLAGS} ${CMAKE_CPP_FLAGS_${CMAKE_BUILD_TYPE}}")

SET(NC_FLIBS "")

get_target_property(ALL_TLL_LIBS netcdff LINK_LIBRARIES)

FOREACH(_LIB IN LISTS ALL_TLL_LIBS)
  GET_FILENAME_COMPONENT(_LIB_NAME ${_LIB} NAME_WE)
  STRING(REGEX REPLACE "^lib" "" _NAME ${_LIB_NAME})
  LIST(APPEND NC_FLIBS "-l${_NAME}")
  GET_FILENAME_COMPONENT(_LIB_DIR ${_LIB} PATH)
  LIST(APPEND LINKFLAGS "-L${_LIB_DIR}")
ENDFOREACH()

###
# Perform some cleanup for the nf-config related variables.
###
IF(NC_FLIBS)
  STRING(REPLACE ";" " " NC_FLIBS "${NC_FLIBS}")
  STRING(REPLACE "-lhdf5::hdf5-shared" "-lhdf5" NC_FLIBS ${NC_FLIBS})
  STRING(REPLACE "-lhdf5::hdf5_hl-shared" "-lhdf5_hl" NC_FLIBS ${NC_FLIBS})
  STRING(REPLACE "-lhdf5::hdf5-static" "-lhdf5" NC_FLIBS ${NC_FLIBS})
  STRING(REPLACE "-lhdf5::hdf5_hl-static" "-lhdf5_hl" NC_FLIBS ${NC_FLIBS})
  STRING(REPLACE "-lnetCDF::netcdf" "-lnetcdf" NC_FLIBS ${NC_FLIBS})
  STRING(REPLACE "-lMPI::MPI_Fortran" "" NC_FLIBS ${NC_FLIBS})
  STRING(REPLACE "-lnetcdff_c" "" NC_FLIBS ${NC_FLIBS})
ENDIF()

STRING(REPLACE ";" " " LINKFLAGS "${LINKFLAGS}")

LIST(REMOVE_DUPLICATES NC_FLIBS)
LIST(REMOVE_DUPLICATES LINKFLAGS)

SET(LIBS ${NC_FLIBS})
SET(NC_FLIBS "-lnetcdff ${NC_FLIBS}")


is_disabled(BUILD_SHARED_LIBS enable_static)
is_enabled(BUILD_SHARED_LIBS enable_shared)

is_enabled(BUILD_F03 HAS_F03)
is_enabled(ENABLE_DAP HAS_DAP)
is_enabled(USE_NETCDF_V2 HAS_NC2)
is_enabled(USE_NETCDF4 HAS_NC4)
is_enabled(HAVE_SZIP_WRITE HAS_SZIP_WRITE)
is_enabled(NF_HAVE_ZSTD HAS_ZSTD)
is_enabled(HAVE_QUANTIZE HAS_QUANTIZE)
is_enabled(USE_LOGGING HAS_LOGGING)
is_enabled(ENABLE_CDF5 HAS_CDF5)
is_enabled(BUILD_PARALLEL HAS_PARALLEL)
is_enabled(NF_HAS_PARALLEL4 HAS_NC4_PARALLEL)
is_enabled(NF_HAS_PNETCDF HAS_PNETCDF)
# Generate file from template.
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/libnetcdff.settings.in"
  "${CMAKE_CURRENT_BINARY_DIR}/libnetcdff.settings"
  @ONLY)

# Read in settings file, print out.
# Avoid using system-specific calls so that this
# might also work on Windows.
FILE(READ "${CMAKE_CURRENT_BINARY_DIR}/libnetcdff.settings"
  LIBNETCDFF_SETTINGS)
MESSAGE(${LIBNETCDFF_SETTINGS})
IF(NC_HAVE_ZSTD)
  IF(NOT NF_HAVE_ZSTD)
    MESSAGE(STATUS "ATTENTION: NO HDF5_PLUGIN_PATH SET")
    IF(NOT DEFINED ENV{HDF5_PLUGIN_PATH})
      MESSAGE(WARNING "libnetcdf was built with zstd support, but HDF5_PLUGIN_PATH is not set.  Please set the HDF5_PLUGIN_PATH environmental variable if zstandard tests fail.")
      MESSAGE(WARNING "You can use `nc-config --plugindir` to determine the appropriate value for HDF5_PLUGIN_PATH")
    ENDIF()
  ENDIF(NOT NF_HAVE_ZSTD)
ENDIF(NC_HAVE_ZSTD)


# Install libnetcdf.settings file into same location
# as the libraries.
INSTALL(FILES "${CMAKE_BINARY_DIR}/libnetcdff.settings"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  COMPONENT libraries)

#####
# End libnetcdff.settings section.
#####

# Set variables used in nf-config and pkgconfig templates
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}")
set(libdir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(includedir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
set(fmoddir "${CMAKE_INSTALL_PREFIX}/${Fortran_INSTALL_MODDIR}")
set(CC "${CMAKE_C_COMPILER}")
set(FC "${CMAKE_Fortran_COMPILER}")
set(CPPFLAGS "${CMAKE_C_FLAGS} ${CMAKE_CPP_FLAGS}")

# Generate nf-config from template
configure_file("${CMAKE_SOURCE_DIR}/nf-config.in"
  "${CMAKE_BINARY_DIR}/tmp/nf-config" @ONLY
  NEWLINE_STYLE LF)
file(COPY "${CMAKE_BINARY_DIR}/tmp/nf-config"
  DESTINATION ${CMAKE_BINARY_DIR}/
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

configure_file(
        ${CMAKE_SOURCE_DIR}/netcdf-fortran.pc.in
        ${CMAKE_BINARY_DIR}/netcdf-fortran.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/netcdf-fortran.pc
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
              COMPONENT utilities)

install(PROGRAMS ${NC4F_BINARY_DIR}/nf-config
        DESTINATION "${CMAKE_INSTALL_BINDIR}"
        COMPONENT utilities)

include(CMakePackageConfigHelpers)

set(ConfigPackageLocation "${CMAKE_INSTALL_LIBDIR}/cmake/netCDF-Fortran")

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/netCDF-FortranConfigVersion.cmake"
  VERSION ${NC4F_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install(
  EXPORT netcdffTargets
  FILE netcdffTargets.cmake
  NAMESPACE netCDF::
  DESTINATION ${ConfigPackageLocation}
  )

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/netCDF-FortranConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/netCDF-FortranConfig.cmake"
  INSTALL_DESTINATION "${ConfigPackageLocation}"
  )


install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/netCDF-FortranConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/netCDF-FortranConfigVersion.cmake"
  DESTINATION "${ConfigPackageLocation}"
  COMPONENT headers
  )
