module AdvSchemeEnumModule
  use KindModule, only: I4B

  implicit none

  ! Advection scheme codes
  integer(I4B), parameter :: ADV_SCHEME_UPSTREAM = 0
  integer(I4B), parameter :: ADV_SCHEME_CENTRAL = 1
  integer(I4B), parameter :: ADV_SCHEME_TVD = 2
  integer(I4B), parameter :: ADV_SCHEME_UTVD = 3
  integer(I4B), parameter :: ADV_SCHEME_ULTIMATE = 4

end module AdvSchemeEnumModule
