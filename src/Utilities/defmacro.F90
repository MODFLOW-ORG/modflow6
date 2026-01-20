module DefinedMacros
  ! -- modules
  use KindModule, only: LGP, I4B
  use ConstantsModule, only: OSUNDEF, OSLINUX, OSMAC, OSWIN
  implicit none
  private
  public :: get_os, is_extended, using_petsc, using_netcdf
contains

  !> @brief Get operating system
  !!
  !! Function to get the integer operating system enum value.
  !!
  !! @return      ios                operating system enum
  !<
  function get_os() result(ios)
    ! -- local variables
    integer(I4B) :: ios !< operating system
    !
    ! -- initialize ios
    ios = OSUNDEF
    !
    ! -- set operating system variables
#ifdef __GFORTRAN__
# ifdef __linux__
    ios = OSLINUX
# endif
# ifdef __APPLE__
    ios = OSMAC
# endif
# ifdef _WIN32
    ios = OSWIN
# endif
#endif
#ifdef __INTEL_COMPILER
# ifdef __linux__
    ios = OSLINUX
# endif
# ifdef __APPLE__
    ios = OSMAC
# endif
# ifdef _WIN32
    ios = OSWIN
# endif
#endif
    !
  end function get_os

  !> @brief Determine if this is the extended version
  !!
  !! Function to get a logical indicating if this is the
  !! extended version of MODFLOW.
  !!
  !! @return      isextended           extended version logical
  !<
  function is_extended() result(isextended)
    ! -- return variables
    logical(LGP) :: isextended !< extended version logical
    ! -- local variables
    logical(LGP) :: ispetsc
    logical(LGP) :: isnetcdf
    !
    ! -- initialize isextended
    isextended = .FALSE.
    !
    ! -- check if using petsc
    ispetsc = using_petsc()
    !
    ! -- check if using netcf
    isnetcdf = using_netcdf()
    !
    !
    if (ispetsc .EQV. .TRUE. .OR. isnetcdf .EQV. .TRUE.) then
      isextended = .TRUE.
    end if
    !
  end function is_extended

  !> @brief Determine if using petsc
  !!
  !! Function to get a logical indicating if petsc is
  !! being used.
  !!
  !! @return      petscused                petsc used logical
  !<
  function using_petsc() result(petscused)
    ! -- return variable
    logical(LGP) :: petscused !< petsc used logical
    !
    ! -- initialize petscavail
    petscused = .FALSE.
    !
    ! -- set operating system variables
#ifdef __WITH_PETSC__
    petscused = .TRUE.
#endif
    !
  end function using_petsc

  !> @brief Determine if using netcdf
  !!
  !! Function to get a logical indicating if netcdf is
  !! being used.
  !!
  !! @return      netcdfused                netcdf used logical
  !<
  function using_netcdf() result(netcdfused)
    ! -- return variable
    logical(LGP) :: netcdfused !< netcdf used logical
    !
    ! -- initialize petscavail
    netcdfused = .FALSE.
    !
    ! -- set operating system variables
#ifdef __WITH_NETCDF__
    netcdfused = .TRUE.
#endif
    !
  end function using_netcdf

end module DefinedMacros
