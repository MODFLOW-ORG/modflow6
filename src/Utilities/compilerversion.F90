module CompilerVersion
  ! -- modules
  use iso_fortran_env, only: compiler_options, compiler_version
  use ConstantsModule, only: LENBIGLINE, &
                             CUNKNOWN, CGFORTRAN, CINTEL, CCRAYFTN
  use KindModule, only: I4B
  implicit none
  private
  ! -- compiler version
  character(len=10) :: c_compiler !< compiler string
  character(len=10) :: c_version !< compiler version string
  character(len=20) :: c_date !< compilation date
  integer(I4B) :: icompiler = CUNKNOWN !< compiler enum
  public :: get_compiler, get_compile_date, get_compile_options
contains

  !> @ brief Get compiler information
  !!
  !!  Subroutine returns a string with compilation date and compiler.
  !!
  !<
  subroutine get_compiler(txt)
    ! -- dummy variables
    character(len=LENBIGLINE), intent(inout) :: txt !< compiler information
    !
    ! -- set variables
#ifdef __GFORTRAN__
    icompiler = CGFORTRAN
    c_date = __DATE__//' '//__TIME__
#endif
#ifdef __INTEL_COMPILER
    icompiler = CINTEL
    c_date = __DATE__//' '//__TIME__
#endif
#ifdef _CRAYFTN
    icompiler = CCRAYFTN
    c_date = __DATE__//' '//__TIME__
#endif
    !
    ! -- set compiler strings
    if (icompiler == CUNKNOWN) then
      c_compiler = 'UNKNOWN'
      c_version = '??.??'
      c_date = '??? ?? ???? ??:??:??'
    end if
    !
    ! -- write string with compiler information
    write (txt, '(a,3(1x,a))') &
      'MODFLOW 6 compiled', trim(adjustl(c_date)), &
      'with', trim(adjustl(compiler_version()))
  end subroutine get_compiler

  !> @ brief Get compilation date
  !!
  !!  Subroutine returns a string with compilation date
  !!
  !<
  subroutine get_compile_date(txt)
    ! -- dummy variables
    character(len=20), intent(inout) :: txt !< compilation date
    ! -- set variables
#ifdef __GFORTRAN__
    c_date = __DATE__//' '//__TIME__
#endif
#ifdef __INTEL_COMPILER
    c_date = __DATE__//' '//__TIME__
#endif
#ifdef _CRAYFTN
    c_date = __DATE__//' '//__TIME__
#endif
    !
    ! -- write compilation date string
    write (txt, '(a)') trim(adjustl(c_date))
  end subroutine get_compile_date

  !> @ brief Get compilation options
  !!
  !!  Subroutine returns a string with compilation options
  !!
  !<
  subroutine get_compile_options(txt)
    ! -- dummy variables
    character(len=LENBIGLINE), intent(inout) :: txt !< compilation options
    ! -- set variables
    !
    ! -- set txt string
    write (txt, '(a)') &
      'MODFLOW 6 compiler options:'//' '//trim(adjustl(compiler_options()))
  end subroutine get_compile_options

end module CompilerVersion
