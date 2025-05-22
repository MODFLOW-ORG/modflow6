!> @brief This module contains version information
!!
!! This module contains subroutines for accessing version information
!! and writing header information to the program listing files.
!!
!<
module VersionModule
  ! -- module imports
  use KindModule
  use DefinedMacros, only: is_extended, using_petsc, using_netcdf
  use ConstantsModule, only: LENBIGLINE, LENHUGELINE, DZERO
  use SimVariablesModule, only: istdout
  use MessageModule, only: write_message, write_message_centered, &
                           write_message_counter
  use CompilerVersion, only: get_compiler, get_compile_options
  implicit none
  public
  ! -- modflow 6 version
  integer(I4B), parameter :: IDEVELOPMODE = 1
  character(len=*), parameter :: VERSIONNUMBER = '6.7.0.dev2'
  character(len=*), parameter :: VERSIONTAG = ' (preliminary) 05/12/2025'
  character(len=40), parameter :: VERSION = VERSIONNUMBER//VERSIONTAG
  character(len=2), parameter :: MFVNAM = ' 6'
  character(len=*), parameter :: MFTITLE = &
    &'U.S. GEOLOGICAL SURVEY MODULAR HYDROLOGIC MODEL'
  character(len=*), parameter :: FMTLICENSE = &
    "(/,&
    &'As a work of the United States Government, this USGS product is ',/,&
    &'in the public domain within the United States. You can copy, ',/,&
    &'modify, distribute, and perform the work, even for commercial ',/,&
    &'purposes, all without asking permission. Additionally, USGS ',/,&
    &'waives copyright and related rights in the work worldwide ',/,&
    &'through CC0 1.0 Universal Public Domain Dedication ',/,&
    &'(https://creativecommons.org/publicdomain/zero/1.0/).',//,&
    &'The following GNU Lesser General Public License (LGPL) libraries',/,&
    &'are used in this USGS product:',//,"// &
    "'    SPARSKIT version 2.0',/,&
    &'      ilut, luson, and qsplit ',/,&
    &'      (https://www-users.cse.umn.edu/~saad/software/SPARSKIT/)',//,&
    &'    RCM - Reverse Cuthill McKee Ordering',/,&
    &'      (https://people.math.sc.edu/Burkardt/f_src/rcm/rcm.html)',//,&
    &'    BLAS - Basic Linear Algebra Subprograms Level 1',/,&
    &'      (https://people.math.sc.edu/Burkardt/f_src/blas1_d/',&
    &'blas1_d.html)',//,"// &
    "'    SPARSEKIT - Sparse Matrix Utility Package',/,&
    &'      amux, dperm, dvperm, rperm, and cperm',/,&
    &'      (https://people.sc.fsu.edu/~jburkardt/f77_src/sparsekit/',&
    &'sparsekit.html)',//,&
    &'The following BSD-3 License libraries are used in this USGS product:',//,&
    &'    Modern Fortran DAG Library',/,&
    &'      Copyright (c) 2018, Jacob Williams',/,&
    &'      All rights reserved.',/,&
    &'      (https://github.com/jacobwilliams/daglib)',/&
    &)"
  character(len=*), parameter :: PETSCLICENSE = &
  "(&
  &'The following 2-clause BSD License library is used in this',/,&
  &'USGS product:',//,&
  &'    PETSc, the Portable, Extensible Toolkit for Scientific',/,&
  &'    Computation Library',/,&
  &'      Copyright (c) 1991-2021, UChicago Argonne, LLC',/,&
  &'      and the PETSc Development Team All rights reserved.',/,&
  &'      (https://petsc.org/release/)',/&
  &)"
  character(len=*), parameter :: NETCDFLICENSE = &
  "(&
  &'The following library is used in this USGS product:',//,&
  &'    NetCDF, network Common Data Form software library',/,&
  &'      Copyright (c) 1993-2014 University Corporation for Atmospheric',/,&
  &'      Research/Unidata. Redistribution and use in source and binary',/,&
  &'      forms, with or without modification, are permitted provided that',/,&
  &'      the conditions in the NetCDF copyright are met',/,&
  &'        (https://www.unidata.ucar.edu/software/netcdf/copyright.html)',/&
  &)"
  ! -- disclaimer must be appropriate for version (release or release candidate)
  character(len=*), parameter :: FMTDISCLAIMER = &
    "(/,&
    &'This software is preliminary or provisional and is subject to ',/,&
    &'revision. It is being provided to meet the need for timely best ',/,&
    &'science. The software has not received final approval by the U.S. ',/,&
    &'Geological Survey (USGS). No warranty, expressed or implied, is made ',/,&
    &'by the USGS or the U.S. Government as to the functionality of the ',/,&
    &'software and related material nor shall the fact of release ',/,&
    &'constitute any such warranty. The software is provided on the ',/,&
    &'condition that neither the USGS nor the U.S. Government shall be held ',/,&
    &'liable for any damages resulting from the authorized or unauthorized ',/,&
    &'use of the software.',/)"

contains

  !> @ brief Write program header
  !!
  !!  Write header for program to the program listing file.
  !!
  !<
  subroutine write_listfile_header(iout, cmodel_type, write_sys_command, &
                                   write_kind_info)
    ! -- dummy variables
    integer(I4B), intent(in) :: iout !< program listing file
    character(len=*), intent(in), optional :: cmodel_type !< optional model type string
    logical(LGP), intent(in), optional :: write_sys_command !< boolean indicating if the system command should be written
    logical(LGP), intent(in), optional :: write_kind_info !< boolean indicating in program data types should be written
    ! -- local variables
    integer(I4B), parameter :: iheader_width = 80
    character(len=22) :: cheader
    character(len=LENBIGLINE) :: syscmd
    character(len=LENBIGLINE) :: compiler
    character(len=LENBIGLINE) :: compiler_options
    logical(LGP) :: wki
    logical(LGP) :: wsc
    !
    ! -- set pro string
    if (is_extended()) then
      write (cheader, '(3a)') 'MODFLOW', MFVNAM, ' EXTENDED'
    else
      write (cheader, '(2a)') 'MODFLOW', MFVNAM
    end if
    !
    ! -- Write title to iout
    call write_message_centered(text=cheader, linelen=iheader_width, &
                                iunit=iout)
    call write_message_centered(text=MFTITLE, linelen=iheader_width, &
                                iunit=iout)
    !
    ! -- Write model type to list file
    if (present(cmodel_type)) then
      call write_message_centered(text=cmodel_type, linelen=iheader_width, &
                                  iunit=iout)
    end if
    !
    ! -- Write version
    call write_message_centered(text='VERSION '//VERSION, &
                                linelen=iheader_width, iunit=iout)
    !
    ! -- Write if develop mode
    if (IDEVELOPMODE == 1) then
      call write_message_centered(text='***DEVELOP MODE***', &
                                  linelen=iheader_width, iunit=iout)
    end if
    !
    ! -- Write compiler version
    call get_compiler(compiler)
    call write_message_centered(text=' ', linelen=iheader_width, &
                                iunit=iout)
    call write_message_centered(text=trim(adjustl(compiler)), &
                                linelen=iheader_width, iunit=iout)
    !
    ! -- Write disclaimer
    write (iout, FMTDISCLAIMER)
    !
    ! -- Write license information
    if (iout /= istdout) then
      call write_license(iout)
    end if
    !
    ! -- write compiler options
    if (iout /= istdout) then
      call get_compile_options(compiler_options)
      call write_message_counter(text=compiler_options, iunit=iout)
    end if
    !
    ! -- Write the system command used to initiate simulation
    wsc = .true.
    if (present(write_sys_command)) wsc = write_sys_command
    if (wsc) then
      call GET_COMMAND(syscmd)
      write (iout, '(/,a,/,a)') &
        'System command used to initiate simulation:', trim(syscmd)
    end if
    !
    ! -- Write precision of real variables
    wki = .true.
    if (present(write_kind_info)) wki = write_kind_info
    if (wki) then
      write (iout, '(/,a)') 'MODFLOW was compiled using uniform precision.'
      call write_kindinfo(iout)
    end if
    write (iout, *)

  end subroutine write_listfile_header

  !> @ brief Write program license
  !!
  !!  Write license for program to the program listing file.
  !!
  !<
  subroutine write_license(iout)
    ! -- dummy variables
    integer(I4B), intent(in), optional :: iout !< program listing file
    !
    ! - write standard license
    if (present(iout)) then
      write (iout, FMTLICENSE)
    else
      call write_message('', fmt=FMTLICENSE)
    end if
    !
    ! -- write NetCDF license
    if (using_netcdf()) then
      if (present(iout)) then
        write (iout, NETCDFLICENSE)
      else
        call write_message('', fmt=NETCDFLICENSE)
      end if
    end if
    !
    ! -- write PETSc license
    if (using_petsc()) then
      if (present(iout)) then
        write (iout, PETSCLICENSE)
      else
        call write_message('', fmt=PETSCLICENSE)
      end if
    end if

  end subroutine write_license

end module VersionModule

