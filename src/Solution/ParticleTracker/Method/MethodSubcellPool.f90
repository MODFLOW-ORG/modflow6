!> @brief Subcell-level tracking methods.
module MethodSubcellPoolModule

  use MethodSubcellPollockModule
  use MethodSubcellTernaryModule
  implicit none

  private
  public :: create_method_subcell_pool
  public :: destroy_method_subcell_pool

  type(MethodSubcellPollockType), pointer, public :: method_subcell_plck => null()
  type(MethodSubcellTernaryType), pointer, public :: method_subcell_tern => null()

contains

  !> @brief Create the subcell method pool
  subroutine create_method_subcell_pool()
    call create_method_subcell_pollock(method_subcell_plck)
    call create_method_subcell_ternary(method_subcell_tern)
  end subroutine create_method_subcell_pool

  !> @brief Destroy the subcell method pool
  subroutine destroy_method_subcell_pool()
    call method_subcell_plck%deallocate()
    deallocate (method_subcell_plck)
    call method_subcell_tern%deallocate()
    deallocate (method_subcell_tern)
  end subroutine destroy_method_subcell_pool

end module MethodSubcellPoolModule
