module SimVariablesPHMFModule
  use ConstantsModule, only: MAXCHARLEN
  use SimPHMFModule, only: store_error, ustop
  implicit none
  
  private
  public :: echo, get_ntimes_global, set_ntimes_global

  logical :: echo = .false.
  integer :: ntimesglobal = -1

contains

  subroutine set_ntimes_global(ntimes)
    ! dummy
    integer, intent(in) :: ntimes
    ! local
    character(len=MAXCHARLEN) :: ermsg
    !
    if (ntimesglobal == -1) then
      ntimesglobal = ntimes
    elseif (ntimesglobal == ntimes) then
      continue
    else
      ermsg = 'Number of observation times is inconsistent.'
      call store_error(ermsg)
      call ustop()
    endif
    ! 
    return
  end subroutine set_ntimes_global
  
  integer function get_ntimes_global()
    if (ntimesglobal == -1) then
      get_ntimes_global = 1
    else
      get_ntimes_global = ntimesglobal
    endif
    return
  end function get_ntimes_global

end module SimVariablesPHMFModule
