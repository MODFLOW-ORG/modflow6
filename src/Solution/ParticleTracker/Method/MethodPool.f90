!> @brief Model-level tracking methods.
module MethodPoolModule
  use MethodModule, only: MethodType
  use MethodDisModule, only: MethodDisType, create_method_dis
  use MethodDisvModule, only: MethodDisvType, create_method_disv
  implicit none

  private
  public :: create_method_pool, destroy_method_pool

  type(MethodDisType), pointer, public :: method_dis => null()
  type(MethodDisvType), pointer, public :: method_disv => null()

contains

  !> @brief Create the method pool
  subroutine create_method_pool()
    call create_method_dis(method_dis)
    call create_method_disv(method_disv)
  end subroutine create_method_pool

  !> @brief Destroy the method pool
  subroutine destroy_method_pool()
    call method_dis%deallocate()
    deallocate (method_dis)
    call method_disv%deallocate()
    deallocate (method_disv)
  end subroutine destroy_method_pool

end module MethodPoolModule
