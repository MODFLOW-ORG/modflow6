!> @brief Cell-level tracking methods.
module MethodCellPoolModule

  use MethodCellPollockModule
  use MethodCellPollockQuadModule
  use MethodCellTernaryModule
  use MethodCellPassToBotModule
  implicit none

  private
  public :: create_method_cell_pool
  public :: destroy_method_cell_pool

  type(MethodCellPollockType), pointer, public :: method_cell_plck => null()
  type(MethodCellPollockQuadType), pointer, public :: method_cell_quad => null()
  type(MethodCellTernaryType), pointer, public :: method_cell_tern => null()
  type(MethodCellPassToBotType), pointer, public :: method_cell_ptb => null()

contains

  !> @brief Create the cell method pool
  subroutine create_method_cell_pool()
    call create_method_cell_pollock(method_cell_plck)
    call create_method_cell_quad(method_cell_quad)
    call create_method_cell_ternary(method_cell_tern)
    call create_method_cell_ptb(method_cell_ptb)
  end subroutine create_method_cell_pool

  !> @brief Destroy the cell method pool
  subroutine destroy_method_cell_pool()
    call method_cell_plck%deallocate()
    deallocate (method_cell_plck)
    call method_cell_quad%deallocate()
    deallocate (method_cell_quad)
    call method_cell_tern%deallocate()
    deallocate (method_cell_tern)
    call method_cell_ptb%deallocate()
    deallocate (method_cell_ptb)
  end subroutine destroy_method_cell_pool

end module MethodCellPoolModule
