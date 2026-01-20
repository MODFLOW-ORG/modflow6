module CharacterStringModule

  implicit none
  private
  public :: CharacterStringType

  !> This class is used to store a single deferred-length
  !! character string.  It was designed to work in an
  !! array implementation so that a jagged character array
  !! could be used in MODFLOW and stored in the memory
  !! manager.
  !!
  !! The overloaded methods allow instances to behave like
  !! a regular string and work with intrinsic Fortran
  !! character strings.  Ideas for the implementation were
  !! inspired by:
  !!   https://gitlab.com/everythingfunctional/iso_varying_string
  !
  !! Can be improved as necessary to overload other string
  !! functions, such as write_formatted, trim, len, ...
  !!
  !<
  type :: CharacterStringType
    private
    character(len=:), allocatable :: charstring
  contains
    procedure, pass(lhs) :: assign_to_charstring
    procedure, pass(rhs) :: assign_from_charstring
    procedure, pass(rhs) :: character_eq_charstring
    procedure, pass(lhs) :: charstring_eq_character
    procedure :: charstring_eq_charstring
    procedure :: write_unformatted
    procedure :: strlen
    procedure :: destroy
    generic :: assignment(=) => assign_to_charstring, assign_from_charstring
    generic :: operator(==) => character_eq_charstring, &
      charstring_eq_character, &
      charstring_eq_charstring
    ! not supported by gfortran 5 and 6
    ! disable for now
    ! generic :: write (unformatted) => write_unformatted
  end type CharacterStringType

contains

  recursive subroutine assign_to_charstring(lhs, rhs)
    class(CharacterStringType), intent(out) :: lhs
    character(len=*), intent(in) :: rhs
    logical :: allocate_charstring
    allocate_charstring = .false.
    if (allocated(lhs%charstring)) then
      if (len(lhs%charstring) <= len(rhs)) then
        lhs%charstring(:) = rhs
      else
        allocate_charstring = .true.
      end if
    else
      allocate_charstring = .true.
    end if
    if (allocate_charstring) then
      lhs%charstring = rhs
    end if
  end subroutine assign_to_charstring

  subroutine assign_from_charstring(lhs, rhs)
    character(len=*), intent(out) :: lhs
    class(CharacterStringType), intent(in) :: rhs
    if (allocated(rhs%charstring)) then
      lhs = rhs%charstring
    else
      lhs = ''
    end if
  end subroutine assign_from_charstring

  elemental function character_eq_charstring(lhs, rhs) result(equals)
    character(len=*), intent(in) :: lhs
    class(CharacterStringType), intent(in) :: rhs
    logical :: equals
    if (allocated(rhs%charstring)) then
      equals = lhs == rhs%charstring
    else
      equals = lhs == ''
    end if
  end function character_eq_charstring

  elemental function charstring_eq_character(lhs, rhs) result(equals)
    class(CharacterStringType), intent(in) :: lhs
    character(len=*), intent(in) :: rhs
    logical :: equals
    if (allocated(lhs%charstring)) then
      equals = lhs%charstring == rhs
    else
      equals = rhs == ''
    end if
  end function charstring_eq_character

  elemental function charstring_eq_charstring(this, rhs) result(equals)
    class(CharacterStringType), intent(in) :: this
    class(CharacterStringType), intent(in) :: rhs
    logical :: equals

    equals = .false.
    if (allocated(this%charstring)) then
      equals = (rhs == this%charstring)
    end if

  end function charstring_eq_charstring

  subroutine write_unformatted(this, unit, iostat, iomsg)
    class(CharacterStringType), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg
    iostat = 0
    if (allocated(this%charstring)) then
      write (unit, iostat=iostat) this%charstring
    end if
  end subroutine write_unformatted

  function strlen(this) result(length)
    class(CharacterStringType), intent(in) :: this
    integer :: length

    if (allocated(this%charstring)) then
      length = len(this%charstring)
    else
      length = 0
    end if
  end function strlen

  subroutine destroy(this)
    class(CharacterStringType), intent(inout) :: this
    if (allocated(this%charstring)) deallocate (this%charstring)
  end subroutine destroy

end module CharacterStringModule
