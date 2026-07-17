module PetscImsPreconditionerModule
#include <petsc/finclude/petscksp.h>
  use petscksp
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMEMPATH, LINELENGTH, DZERO, DONE, IZERO
  use PetscMatrixModule, only: PetscMatrixType
  use SparseMatrixModule, only: SparseMatrixType
  use ImsLinearSettingsModule, only: ImsLinearSettingsType, &
                                     IPC_ILU0, IPC_MILU0, IPC_ILUT, &
                                     IPC_MILUT, IPC_AMG
  use ImsLinearAmgModule, only: ImsAmgDataType, ims_amg_setup, &
                                ims_amg_apply, ims_amg_da
  use SimModule, only: store_error
  use SimVariablesModule, only: iout

  implicit none
  private

  public :: PcShellCtxType
  public :: pcshell_apply
  public :: pcshell_setup

  ! preconditioner name indexed by the IPC_* type (IPC_ILU0=1 ... IPC_AMG=5)
  character(len=9), dimension(5) :: IMS_PC_TYPE = ["IMS ILU0 ", &
                                                   "IMS MILU0", &
                                                   "IMS ILUT ", &
                                                   "IMS MILUT", &
                                                   "IMS AMG  "]

  type :: PcShellCtxType
    character(len=LENMEMPATH) :: memory_path !< the memory path, taken from the system matrix
    class(PetscMatrixType), pointer :: system_matrix !< pointer to the petsc system matrix
    type(ImsLinearSettingsType), pointer :: linear_settings !< the linear settings from IMS
    integer(I4B) :: ipc !< the preconditioner type (an IPC_* enum: IPC_ILU0, IPC_MILU0, IPC_ILUT, IPC_MILUT, or IPC_AMG)
    character(len=16) :: ims_pc_type !< the preconditioner type as a string
    ! PC matrix
    integer(I4B), dimension(:), contiguous, pointer :: IAPC => null() !< CSR row pointer
    integer(I4B), dimension(:), contiguous, pointer :: JAPC => null() !< CSR columns
    real(DP), dimension(:), contiguous, pointer :: APC => null() !< CSR Preconditioner matrix values
    ! work vectors
    integer(I4B), dimension(:), contiguous, pointer :: IW => null() !< (M)ILU0 work array
    real(DP), dimension(:), contiguous, pointer :: W => null() !< (M)ILU0 work array
    integer(I4B), dimension(:), contiguous, pointer :: JLU => null() !< (M)ILUT work array
    integer(I4B), dimension(:), contiguous, pointer :: JW => null() !< (M)ILUT work array
    real(DP), dimension(:), contiguous, pointer :: WLU => null() !< (M)ILUT work array
    ! AMG hierarchy (block-local, allocated when ipc == IPC_AMG)
    type(ImsAmgDataType), pointer :: amg => null() !< block-local AMG hierarchy
  contains
    procedure :: create => pctx_create
    procedure :: destroy => pctx_destroy

  end type

  interface
    subroutine PCShellGetContext(pc, ctx, ierr)
      import PcShellCtxType, tPC
      type(tPC) :: pc
      type(PcShellCtxType), pointer :: ctx
      integer :: ierr
    end subroutine
  end interface

contains

  !> @brief Create the preconditioner context from the system matrix
  !<
  subroutine pctx_create(this, matrix, lin_settings)
    use IMSLinearBaseModule, only: ims_base_pccrs, ims_calc_pcdims
    use MemoryManagerModule, only: mem_allocate
    class(PcShellCtxType) :: this !< this instance
    class(PetscMatrixType), pointer :: matrix !< the linear system matrix
    type(ImsLinearSettingsType), pointer :: lin_settings !< linear settings from IMS
    ! local
    integer(I4B) :: n, neq, nja
    integer(I4B), dimension(:), contiguous, pointer :: ia, ja
    real(DP), dimension(:), contiguous, pointer :: amat
    integer(I4B) :: niapc, njapc, njlu, njw, nwlu
    character(len=LINELENGTH) :: errmsg

    this%memory_path = matrix%memory_path
    this%linear_settings => lin_settings
    this%system_matrix => matrix

    ! -- determine the preconditioner type; AMG (PRECONDITIONER_TYPE AMG) is
    !    applied per subdomain (block-local), just like the MILU0/MILUT ones
    if (lin_settings%ipc_type == 1) then
      this%ipc = IPC_AMG
    else
      this%ipc = IPC_ILU0
      if (this%linear_settings%level > 0) this%ipc = this%ipc + 2
      if (this%linear_settings%relax > 0.0_DP) this%ipc = this%ipc + 1
    end if
    this%ims_pc_type = IMS_PC_TYPE(this%ipc)

    if (this%ipc == IPC_AMG) then
      ! PETSc uses a Krylov accelerator, so the AMG cycle must be a fixed linear
      ! operator. Adaptive omega (RELAXATION_FACTOR 0) makes the preconditioner
      ! variable and is not supported in parallel; require RELAXATION_FACTOR = 1.
      if (lin_settings%relax /= DONE) then
        write (errmsg, '(a)') &
          'PRECONDITIONER_TYPE AMG in a parallel (PETSc) solution requires '// &
          'RELAXATION_FACTOR = 1; adaptive or damped omega is not supported.'
        call store_error(errmsg)
      end if
      allocate (this%amg)
      ! allocate the ILU work arrays at size zero so destroy stays symmetric
      call mem_allocate(this%IAPC, 0, 'IAPC', this%memory_path)
      call mem_allocate(this%JAPC, 0, 'JAPC', this%memory_path)
      call mem_allocate(this%APC, 0, 'APC', this%memory_path)
      call mem_allocate(this%IW, 0, 'IW', this%memory_path)
      call mem_allocate(this%W, 0, 'W', this%memory_path)
      call mem_allocate(this%JLU, 0, 'JLU', this%memory_path)
      call mem_allocate(this%JW, 0, 'JW', this%memory_path)
      call mem_allocate(this%WLU, 0, 'WLU', this%memory_path)
    else
      call matrix%get_aij_local(ia, ja, amat)

      neq = matrix%nrow
      nja = matrix%nnz_local
      call ims_calc_pcdims(neq, nja, ia, this%linear_settings%level, this%ipc, &
                           niapc, njapc, njlu, njw, nwlu)

      ! preconditioner matrix
      call mem_allocate(this%IAPC, niapc + 1, 'IAPC', this%memory_path)
      call mem_allocate(this%JAPC, njapc, 'JAPC', this%memory_path)
      call mem_allocate(this%APC, njapc, 'APC', this%memory_path)
      if (this%ipc == 1 .or. this%ipc == 2) then
        call ims_base_pccrs(niapc, njapc, ia, ja, this%IAPC, this%JAPC)
      end if
      do n = 1, njapc
        this%APC(n) = DZERO
      end do

      ! (M)ILU0 work arrays
      call mem_allocate(this%IW, niapc, 'IW', this%memory_path)
      call mem_allocate(this%W, niapc, 'W', this%memory_path)
      do n = 1, niapc
        this%IW(n) = IZERO
        this%W(n) = DZERO
      end do

      ! (M)ILUT work arrays
      call mem_allocate(this%JLU, njlu, 'JLU', this%memory_path)
      call mem_allocate(this%JW, njw, 'JW', this%memory_path)
      call mem_allocate(this%WLU, nwlu, 'WLU', this%memory_path)

      do n = 1, njlu
        this%JLU(n) = IZERO
      end do
      do n = 1, njw
        this%JW(n) = IZERO
      end do
      do n = 1, nwlu
        this%WLU(n) = DZERO
      end do

    end if

  end subroutine pctx_create

  !> @brief Cleanup the context
  !<
  subroutine pctx_destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    class(PcShellCtxType) :: this

    ! AMG hierarchy
    if (this%ipc == IPC_AMG) then
      call ims_amg_da(this%amg)
      deallocate (this%amg)
      this%amg => null()
    end if

    ! matrix
    call mem_deallocate(this%IAPC)
    call mem_deallocate(this%JAPC)
    call mem_deallocate(this%APC)

    ! work arrays
    call mem_deallocate(this%IW)
    call mem_deallocate(this%W)
    call mem_deallocate(this%JLU)
    call mem_deallocate(this%JW)
    call mem_deallocate(this%WLU)

    this%system_matrix => null()

  end subroutine pctx_destroy

  !> @brief Apply shell preconditioner
  !< (this is not a type bound procedure)
  subroutine pcshell_apply(pc, x, y, ierr)
    use IMSLinearBaseModule, only: ims_base_ilu0a
    external lusol ! from ilut.f90
    PC :: pc !< the shell preconditioner
    Vec :: x !< the input vector
    Vec :: y !< the output vector
    PetscErrorCode :: ierr !< PETSc error code
    ! local
    type(PcShellCtxType), pointer :: pc_ctx => null()
    real(DP), dimension(:), pointer :: local_x, local_y
    integer(I4B) :: neq, nja

    call PCShellGetContext(pc, pc_ctx, ierr)
    CHKERRQ(ierr)

    call VecGetArrayReadF90(x, local_x, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(y, local_y, ierr)
    CHKERRQ(ierr)

    neq = pc_ctx%system_matrix%nrow
    nja = pc_ctx%system_matrix%nnz_local
    if (pc_ctx%ipc == IPC_AMG) then
      call ims_amg_apply(pc_ctx%amg, neq, local_x, local_y)
    else
      select case (pc_ctx%ipc)
      case (IPC_ILU0, IPC_MILU0)
        call ims_base_ilu0a(nja, neq, pc_ctx%APC, pc_ctx%IAPC, pc_ctx%JAPC, &
                            local_x, local_y)
      case (IPC_ILUT, IPC_MILUT)
        call lusol(neq, local_x, local_y, pc_ctx%APC, pc_ctx%JLU, pc_ctx%IW)
      end select
    end if

    call VecRestoreArrayF90(x, local_x, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(y, local_y, ierr)
    CHKERRQ(ierr)

  end subroutine pcshell_apply

  !> @brief Set up the custom preconditioner
  !< (this is not a type bound procedure)
  subroutine pcshell_setup(pc, ierr)
    use IMSLinearBaseModule, only: ims_base_pcu
    PC :: pc !< the shell preconditioner
    PetscErrorCode :: ierr !< PETSc error code
    ! local
    type(PcShellCtxType), pointer :: pc_ctx => null()
    integer(I4B) :: neq, nja
    integer(I4B) :: niapc, njapc, njlu, njw, nwlu
    integer(I4B), dimension(:), contiguous, pointer :: ia, ja
    real(DP), dimension(:), contiguous, pointer :: amat

    call PCShellGetContext(pc, pc_ctx, ierr)
    CHKERRQ(ierr)

    ! note the two different matrix types here:
    ! bridging between PETSc and IMS
    call pc_ctx%system_matrix%get_aij_local(ia, ja, amat)

    neq = pc_ctx%system_matrix%nrow
    nja = pc_ctx%system_matrix%nnz_local

    ! -- AMG builds (or value-updates) the block-local hierarchy; otherwise the
    !    ILU(0)/ILUT factorization is computed below
    if (pc_ctx%ipc == IPC_AMG) then
      call ims_amg_setup(pc_ctx%amg, neq, nja, ia, ja, amat, &
                         max(1, pc_ctx%linear_settings%level), &
                         pc_ctx%linear_settings%nsmooth, &
                         pc_ctx%linear_settings%relax, &
                         pc_ctx%linear_settings%strength_threshold, &
                         pc_ctx%linear_settings%smoother_type)
    else
      niapc = size(pc_ctx%IAPC) - 1
      njapc = size(pc_ctx%JAPC)
      njlu = size(pc_ctx%JLU)
      njw = size(pc_ctx%JW)
      nwlu = size(pc_ctx%WLU)
      call ims_base_pcu(iout, nja, neq, niapc, njapc, &
                        pc_ctx%ipc, pc_ctx%linear_settings%relax, &
                        amat, ia, ja, &
                        pc_ctx%APC, pc_ctx%IAPC, pc_ctx%JAPC, &
                        pc_ctx%IW, pc_ctx%W, &
                        pc_ctx%linear_settings%level, &
                        pc_ctx%linear_settings%droptol, &
                        njlu, njw, nwlu, &
                        pc_ctx%JLU, pc_ctx%JW, pc_ctx%WLU)

    end if

  end subroutine pcshell_setup

end module PetscImsPreconditionerModule
