! Stub for non-Meson builds (MSVS, make, etc.). The real DfnSpec.f90 is
! generated at build time by the Meson custom_target via:
!   python utils/idmloader/scripts/dfn2f90.py --spec <output>
module DfnSpecModule
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  private
  public :: write_spec

contains

  subroutine write_spec()
    write (output_unit, '(a)') &
      'DFN spec not available: rebuild with Meson to enable this feature)'
  end subroutine write_spec

end module DfnSpecModule
