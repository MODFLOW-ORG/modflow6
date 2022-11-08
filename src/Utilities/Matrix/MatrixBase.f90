module MatrixModule
  use ConstantsModule, only: LENMEMPATH
  use KindModule, only: I4B, DP
  use SparseModule, only: sparsematrix
  implicit none
  private

  type, public, abstract :: MatrixBaseType
    character(len=LENMEMPATH) :: memory_path
  contains
    procedure(create_if), deferred :: create
    procedure(destroy_if), deferred :: destroy

    procedure(get_value_pos_if), deferred :: get_value_pos
    procedure(get_diag_value_if), deferred :: get_diag_value

    procedure(set_diag_value_if), deferred :: set_diag_value
    procedure(set_value_pos_if), deferred :: set_value_pos
    procedure(add_value_pos_if), deferred :: add_value_pos
    procedure(add_diag_value_if), deferred :: add_diag_value
    procedure(zero_entries_if), deferred :: zero_entries
    procedure(zero_row_offdiag_if), deferred :: zero_row_offdiag

    procedure(get_first_col_pos_if), deferred :: get_first_col_pos
    procedure(get_last_col_pos_if), deferred :: get_last_col_pos
    procedure(get_column_if), deferred :: get_column
    procedure(get_position_if), deferred :: get_position
    procedure(get_position_diag_if), deferred :: get_position_diag

  end type MatrixBaseType

  abstract interface
    subroutine create_if(this, sparse, mem_path)
      import MatrixBaseType, sparsematrix
      class(MatrixBaseType) :: this
      type(sparsematrix) :: sparse
      character(len=*) :: mem_path
    end subroutine
    subroutine destroy_if(this)
      import MatrixBaseType
      class(MatrixBaseType) :: this
    end subroutine

    function get_value_pos_if(this, ipos) result(value)
      import MatrixBaseType, I4B, DP
      class(MatrixBaseType) :: this
      integer(I4B) :: ipos
      real(DP) :: value
    end function
    function get_diag_value_if(this, irow) result(diag_value)
      import MatrixBaseType, I4B, DP
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      real(DP) :: diag_value
    end function

    subroutine set_diag_value_if(this, irow, diag_value)
      import MatrixBaseType, I4B, DP
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      real(DP) :: diag_value
    end subroutine
    subroutine set_value_pos_if(this, ipos, value)
      import MatrixBaseType, I4B, DP
      class(MatrixBaseType) :: this
      integer(I4B) :: ipos
      real(DP) :: value
    end subroutine
    subroutine add_value_pos_if(this, ipos, value)
      import MatrixBaseType, I4B, DP
      class(MatrixBaseType) :: this
      integer(I4B) :: ipos
      real(DP) :: value
    end subroutine
    subroutine add_diag_value_if(this, irow, value)
      import MatrixBaseType, I4B, DP
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      real(DP) :: value
    end subroutine
    subroutine zero_entries_if(this)
      import MatrixBaseType
      class(MatrixBaseType) :: this
    end subroutine
    subroutine zero_row_offdiag_if(this, irow)
      import MatrixBaseType, I4B
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
    end subroutine
    function get_first_col_pos_if(this, irow) result(first_col_pos)
      import MatrixBaseType, I4B
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      integer(I4B) :: first_col_pos
    end function
    function get_last_col_pos_if(this, irow) result(last_col_pos)
      import MatrixBaseType, I4B
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      integer(I4B) :: last_col_pos
    end function
    function get_column_if(this, ipos) result(icol)
      import MatrixBaseType, I4B
      class(MatrixBaseType) :: this
      integer(I4B) :: ipos
      integer(I4B) :: icol
    end function

    !> @brief Get position index for this (irow,icol) element
    !! in the matrix for direct access with the other routines
    !< Returns -1 when not found.
    function get_position_if(this, irow, icol) result(ipos)
      import MatrixBaseType, I4B
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      integer(I4B) :: icol
      integer(I4B) :: ipos
    end function
    function get_position_diag_if(this, irow) result(ipos_diag)
      import MatrixBaseType, I4B
      class(MatrixBaseType) :: this
      integer(I4B) :: irow
      integer(I4B) :: ipos_diag
    end function
  end interface

end module MatrixModule
