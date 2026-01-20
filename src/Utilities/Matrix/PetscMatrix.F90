module PetscMatrixModule
#include <petsc/finclude/petscksp.h>
  use petscksp
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO
  use SparseModule, only: sparsematrix
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MatrixBaseModule
  use VectorBaseModule
  use PetscVectorModule
  implicit none
  private

  type, public, extends(MatrixBaseType) :: PetscMatrixType
    Mat :: mat !< the PETSc matrix object, NOTE: update() should be called before using this,
               !! in case the matrix CSR array has changed!!!
    integer(I4B) :: nrow !< number of rows in this portion of the global matrix
    integer(I4B) :: ncol !< number of columns in the matrix
    integer(I4B) :: nnz !< number of nonzeros in the matrix
    integer(I4B) :: nnz_local !< number of nonzeros in the local matrix (the diagonal block tied to this process)
    integer(I4B) :: offset !< global offset for the first row in the matrix
    integer(I4B), dimension(:), pointer, contiguous, private :: ia_local !< IA(CSR) for local block, contains 1-based index values
    integer(I4B), dimension(:), pointer, contiguous, private :: ja_local !< JA(CSR) for local block, contains 1-based index values
    real(DP), dimension(:), pointer, contiguous, private :: amat_local !< A(CSR) for local block
    integer(I4B), dimension(:), pointer, contiguous :: ia_petsc !< IA(CSR) for petsc, contains 0-based index values
    integer(I4B), dimension(:), pointer, contiguous :: ja_petsc !< JA(CSR) for petsc, contains 0-based index values
    real(DP), dimension(:), pointer, contiguous :: amat_petsc !< A(CSR) for petsc
    logical(LGP) :: is_parallel
  contains
    ! override
    procedure :: init => pm_init
    procedure :: destroy => pm_destroy
    procedure :: create_vec_mm => pm_create_vec_mm
    procedure :: create_vec => pm_create_vec
    procedure :: get_value_pos => pm_get_value_pos
    procedure :: get_diag_value => pm_get_diag_value
    procedure :: set_diag_value => pm_set_diag_value
    procedure :: set_value_pos => pm_set_value_pos
    procedure :: add_value_pos => pm_add_value_pos
    procedure :: add_diag_value => pm_add_diag_value
    procedure :: zero_entries => pm_zero_entries
    procedure :: zero_row_offdiag => pm_zero_row_offdiag
    procedure :: get_first_col_pos => pm_get_first_col_pos
    procedure :: get_last_col_pos => pm_get_last_col_pos
    procedure :: get_column => pm_get_column
    procedure :: get_position => pm_get_position
    procedure :: get_position_diag => pm_get_position_diag
    procedure :: get_aij => pm_get_aij
    procedure :: get_row_offset => pm_get_row_offset

    procedure :: multiply => pm_multiply

    ! public
    procedure :: update => pm_update
    procedure :: get_aij_local => pm_get_aij_local

    ! private
    procedure, private :: pm_get_position

  end type PetscMatrixType

contains

  subroutine pm_init(this, sparse, mem_path)
    use SimVariablesModule, only: simulation_mode
    class(PetscMatrixType) :: this
    type(sparsematrix) :: sparse
    character(len=*) :: mem_path
    ! local
    PetscErrorCode :: ierr
    integer(I4B) :: i, ierror
    Mat :: local_mat
    integer(I4B), dimension(:), pointer :: ia_tmp, ja_tmp
    integer(I4B) :: nrows_local
    logical(LGP) :: done

    this%memory_path = mem_path
    this%nrow = sparse%nrow
    this%ncol = sparse%ncol
    this%nnz = sparse%nnz
    this%nnz_local = sparse%nnz - sparse%nnz_od
    this%offset = sparse%offset

    ! allocate the diagonal block of the matrix
    call mem_allocate(this%ia_local, this%nrow + 1, 'IA_LOCAL', this%memory_path)
    call mem_allocate(this%ja_local, this%nnz_local, 'JA_LOCAL', this%memory_path)
    call mem_allocate(this%amat_local, this%nnz_local, 'AMAT_LOCAL', &
                      this%memory_path)
    call mem_allocate(this%ia_petsc, this%nrow + 1, 'IA_PETSC', this%memory_path)
    call mem_allocate(this%ja_petsc, this%nnz, 'JA_PETSC', this%memory_path)
    call mem_allocate(this%amat_petsc, this%nnz, 'AMAT_PETSC', this%memory_path)

    call sparse%sort(with_csr=.true.) !< PETSc has full row sorting, MF6 had diagonal first and then sorted
    call sparse%filliaja(this%ia_petsc, this%ja_petsc, ierror)

    ! go to C indexing for PETSc internals
    do i = 1, this%nrow + 1
      this%ia_petsc(i) = this%ia_petsc(i) - 1
    end do
    do i = 1, this%nnz
      this%ja_petsc(i) = this%ja_petsc(i) - 1
      this%amat_petsc(i) = DZERO
    end do

    ! create PETSc matrix object
    if (simulation_mode == 'PARALLEL') then
      this%is_parallel = .true.
      call MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, &
                                     sparse%nrow, sparse%nrow, &
                                     sparse%ncol, sparse%ncol, &
                                     this%ia_petsc, this%ja_petsc, &
                                     this%amat_petsc, this%mat, &
                                     ierr)
    else
      this%is_parallel = .false.
      call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, &
                                     sparse%nrow, sparse%ncol, &
                                     this%ia_petsc, this%ja_petsc, &
                                     this%amat_petsc, this%mat, &
                                     ierr)
    end if
    CHKERRQ(ierr)

    ! set up local block sparsity structure
    nrows_local = 0
    call MatGetDiagonalBlock(this%mat, local_mat, ierr)
    CHKERRQ(ierr)

    call MatGetRowIJF90(local_mat, 1, PETSC_FALSE, PETSC_FALSE, &
                        nrows_local, ia_tmp, ja_tmp, &
                        done, ierr)
    CHKERRQ(ierr)

    do i = 1, size(ia_tmp)
      this%ia_local(i) = ia_tmp(i)
    end do
    do i = 1, size(ja_tmp) - 1
      this%ja_local(i) = ja_tmp(i)
    end do

    call MatRestoreRowIJF90(local_mat, 1, PETSC_FALSE, PETSC_FALSE, &
                            nrows_local, ia_tmp, ja_tmp, done, ierr)
    CHKERRQ(ierr)

  end subroutine pm_init

  !> @brief Copies the values from the CSR array into
  !< the PETSc matrix object
  subroutine pm_update(this)
    class(PetscMatrixType) :: this
    ! local
    PetscErrorCode :: ierr

    if (this%is_parallel) then
      call MatUpdateMPIAIJWithArrays(this%mat, this%nrow, this%nrow, &
                                     this%ncol, this%ncol, &
                                     this%ia_petsc, this%ja_petsc, &
                                     this%amat_petsc, ierr)
      CHKERRQ(ierr)
    end if

  end subroutine pm_update

  subroutine pm_destroy(this)
    class(PetscMatrixType) :: this
    ! local
    PetscErrorCode :: ierr

    call MatDestroy(this%mat, ierr)
    CHKERRQ(ierr)

    call mem_deallocate(this%ia_local)
    call mem_deallocate(this%ja_local)
    call mem_deallocate(this%amat_local)
    call mem_deallocate(this%ia_petsc)
    call mem_deallocate(this%ja_petsc)
    call mem_deallocate(this%amat_petsc)

  end subroutine pm_destroy

  function pm_create_vec_mm(this, n, name, mem_path) result(vec)
    class(PetscMatrixType) :: this
    integer(I4B) :: n !< the nr. of elements in the vector
    character(len=*) :: name !< the variable name (for access through memory manager)
    character(len=*) :: mem_path !< memory path for storing the underlying memory items
    class(VectorBaseType), pointer :: vec !< the vector object to return
    ! local
    class(PetscVectorType), pointer :: petsc_vec

    allocate (petsc_vec)
    call petsc_vec%create_mm(n, name, mem_path)
    vec => petsc_vec

  end function pm_create_vec_mm

  function pm_create_vec(this, n) result(vec)
    class(PetscMatrixType) :: this
    integer(I4B) :: n !< the nr. of elements in the vector
    class(VectorBaseType), pointer :: vec !< the vector object to return
    ! local
    class(PetscVectorType), pointer :: petsc_vec

    allocate (petsc_vec)
    call petsc_vec%create(n)
    vec => petsc_vec

  end function pm_create_vec

  function pm_get_value_pos(this, ipos) result(value)
    class(PetscMatrixType) :: this
    integer(I4B) :: ipos
    real(DP) :: value

    value = this%amat_petsc(ipos)

  end function pm_get_value_pos

  function pm_get_diag_value(this, irow) result(diag_value)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    real(DP) :: diag_value
    ! local
    integer(I4B) :: idiag

    idiag = this%get_position_diag(irow)
    diag_value = this%amat_petsc(idiag)

  end function pm_get_diag_value

  subroutine pm_set_diag_value(this, irow, diag_value)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    real(DP) :: diag_value
    ! local
    integer(I4B) :: idiag

    idiag = this%get_position_diag(irow)
    this%amat_petsc(idiag) = diag_value

  end subroutine pm_set_diag_value

  subroutine pm_set_value_pos(this, ipos, value)
    class(PetscMatrixType) :: this
    integer(I4B) :: ipos
    real(DP) :: value

    this%amat_petsc(ipos) = value

  end subroutine pm_set_value_pos

  subroutine pm_add_value_pos(this, ipos, value)
    class(PetscMatrixType) :: this
    integer(I4B) :: ipos
    real(DP) :: value

    this%amat_petsc(ipos) = this%amat_petsc(ipos) + value

  end subroutine pm_add_value_pos

  subroutine pm_add_diag_value(this, irow, value)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    real(DP) :: value
    ! local
    integer(I4B) :: idiag

    idiag = this%get_position_diag(irow)
    this%amat_petsc(idiag) = this%amat_petsc(idiag) + value

  end subroutine pm_add_diag_value

  function pm_get_first_col_pos(this, irow) result(first_col_pos)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: first_col_pos
    ! local
    integer(I4B) :: irow_local

    ! convert to local row index
    irow_local = irow - this%offset

    ! includes conversion to Fortran's 1-based
    first_col_pos = this%ia_petsc(irow_local) + 1

  end function pm_get_first_col_pos

  function pm_get_last_col_pos(this, irow) result(last_col_pos)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: last_col_pos
    ! local
    integer(I4B) :: irow_local

    ! convert to local row index
    irow_local = irow - this%offset

    ! includes conversion to Fortran's 1-based
    last_col_pos = this%ia_petsc(irow_local + 1)

  end function pm_get_last_col_pos

  function pm_get_column(this, ipos) result(icol)
    class(PetscMatrixType) :: this
    integer(I4B) :: ipos
    integer(I4B) :: icol

    ! includes conversion to Fortran's 1-based
    icol = this%ja_petsc(ipos) + 1

  end function pm_get_column

  !> @brief Return position index for (irow,icol) element
  !! in the matrix. This index can be used in other
  !! routines for direct access.
  !< Returns -1 when not found.
  function pm_get_position(this, irow, icol) result(ipos)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: icol
    integer(I4B) :: ipos
    ! local
    integer(I4B) :: ipos_f
    integer(I4B) :: irow_local

    ipos = -1

    ! convert to local row index
    irow_local = irow - this%offset

    ! includes conversion to Fortran's 1-based
    do ipos_f = this%ia_petsc(irow_local) + 1, this%ia_petsc(irow_local + 1)
      if (this%ja_petsc(ipos_f) + 1 == icol) then
        ipos = ipos_f
        return
      end if
    end do

  end function pm_get_position

  function pm_get_position_diag(this, irow) result(ipos_diag)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: ipos_diag

    ipos_diag = this%pm_get_position(irow, irow)

  end function pm_get_position_diag

  !> @brief Set all entries in the matrix to zero
  !<
  subroutine pm_zero_entries(this)
    class(PetscMatrixType) :: this
    ! local
    integer(I4B) :: i

    do i = 1, this%nnz
      this%amat_petsc(i) = DZERO
    end do

  end subroutine pm_zero_entries

  !> @brief Set all off-diagonal entries in the matrix to zero
  !<
  subroutine pm_zero_row_offdiag(this, irow)
    class(PetscMatrixType) :: this
    integer(I4B) :: irow
    ! local
    integer(I4B) :: ipos, idiag

    idiag = this%get_position_diag(irow)
    do ipos = this%get_first_col_pos(irow), this%get_last_col_pos(irow)
      if (ipos == idiag) cycle
      this%amat_petsc(ipos) = DZERO
    end do

  end subroutine pm_zero_row_offdiag

  subroutine pm_get_aij(this, ia, ja, amat)
    use SimModule, only: ustop
    class(PetscMatrixType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ia
    integer(I4B), dimension(:), pointer, contiguous :: ja
    real(DP), dimension(:), pointer, contiguous :: amat

    call ustop("pm_get_aij not implemented")

  end subroutine pm_get_aij

  !> @brief Get the ia/ja and amat values for the local
  !< diagonal block of this matrix. It's a copy.
  subroutine pm_get_aij_local(this, ia, ja, amat)
    class(PetscMatrixType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ia
    integer(I4B), dimension(:), pointer, contiguous :: ja
    real(DP), dimension(:), pointer, contiguous :: amat
    ! local
    PetscErrorCode :: ierr
    integer(I4B) :: irow, icol, ipos
    integer(I4B), dimension(1) :: irow_idxs, icol_idxs
    real(DP), dimension(1) :: values
    Mat :: local_mat

    call MatGetDiagonalBlock(this%mat, local_mat, ierr)
    CHKERRQ(ierr)

    do irow = 1, this%nrow
      do ipos = this%ia_local(irow), this%ia_local(irow + 1) - 1
        icol = this%ja_local(ipos)
        irow_idxs(1) = irow - 1
        icol_idxs(1) = icol - 1
        call MatGetValues(local_mat, 1, irow_idxs, 1, icol_idxs, values, ierr)
        CHKERRQ(ierr)
        this%amat_local(ipos) = values(1)
      end do
    end do

    ia => this%ia_local
    ja => this%ja_local
    amat => this%amat_local

  end subroutine pm_get_aij_local

  function pm_get_row_offset(this) result(offset)
    class(PetscMatrixType) :: this
    integer(I4B) :: offset

    offset = this%offset

  end function pm_get_row_offset

  !> @brief Calculates global matrix vector product: y = A * x
  !<
  subroutine pm_multiply(this, vec_x, vec_y)
    class(PetscMatrixType) :: this
    class(VectorBaseType), pointer :: vec_x
    class(VectorBaseType), pointer :: vec_y
    ! local
    PetscErrorCode :: ierr
    class(PetscVectorType), pointer :: x, y

    x => null()
    select type (vec_x)
    class is (PetscVectorType)
      x => vec_x
    end select
    y => null()
    select type (vec_y)
    class is (PetscVectorType)
      y => vec_y
    end select

    ! copy data into petsc object
    call this%update()
    ! and multiply
    call MatMult(this%mat, x%vec_impl, y%vec_impl, ierr)
    CHKERRQ(ierr)

  end subroutine pm_multiply

end module PetscMatrixModule
