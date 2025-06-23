module AdvSchemeEnumModule
  use KindModule, only: I4B

  implicit none

  ! Advection scheme codes
  integer(I4B), parameter :: ADV_SCHEME_UPSTREAM = 0
  integer(I4B), parameter :: ADV_SCHEME_CENTRAL = 1
  integer(I4B), parameter :: ADV_SCHEME_TVD = 2

end module AdvSchemeEnumModule