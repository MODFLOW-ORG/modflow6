!> @brief Disable development features in release mode
module FeatureFlagsModule
  use KindModule, only: I4B
  use VersionModule, only: IDEVELOPMODE
  use SimModule, only: store_error, store_error_unit
  implicit none
  private
  public :: developmode, is_release_mode

contains

  !> @brief True if development features should be disabled
  pure logical function is_release_mode()
    is_release_mode = (IDEVELOPMODE == 0)
  end function is_release_mode

  !> @brief Terminate if in release mode (guard development features)
  !!
  !! Terminate the program with an error if the IDEVELOPMODE flag
  !! is set to 0. This allows developing features on the mainline
  !! while disabling them in release builds. An optional file unit
  !! may be specified to associate the feature with an input file.
  !!
  !<
  subroutine developmode(errmsg, iunit)
    ! -- dummy
    character(len=*), intent(in) :: errmsg
    integer(I4B), intent(in), optional :: iunit

    ! -- store error and terminate if in release mode
    if (is_release_mode()) then
      if (present(iunit)) then
        call store_error(errmsg, terminate=.false.)
        call store_error_unit(iunit, terminate=.true.)
      else
        call store_error(errmsg, terminate=.true.)
      end if
    end if

  end subroutine developmode

end module FeatureFlagsModule
