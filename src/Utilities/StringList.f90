module StringListModule

  use KindModule, only: DP, I4B
  use ListModule, only: ListType
  use CharacterStringModule, only: CharacterStringType

  private
  public :: AddStringToList, GetStringFromList

contains

  subroutine ConstructCharacterContainer(newCharCont, text)
    implicit none
    type(CharacterStringType), pointer, intent(out) :: newCharCont
    character(len=*), intent(in) :: text
    !
    allocate (newCharCont)
    newCharCont = text
  end subroutine ConstructCharacterContainer

  function CastAsCharacterStringType(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    type(CharacterStringType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (CharacterStringType)
      res => obj
    end select
  end function CastAsCharacterStringType

  subroutine AddStringToList(list, string)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    character(len=*), intent(in) :: string
    ! -- local
    class(*), pointer :: obj
    type(CharacterStringType), pointer :: newCharacterContainer
    !
    newCharacterContainer => null()
    call ConstructCharacterContainer(newCharacterContainer, string)
    if (associated(newCharacterContainer)) then
      obj => newCharacterContainer
      call list%Add(obj)
    end if
  end subroutine AddStringToList

  function GetStringFromList(list, indx) result(string)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: indx
    character(len=:), allocatable :: string
    ! -- local
    class(*), pointer :: obj
    type(CharacterStringType), pointer :: charcont
    !
    obj => list%GetItem(indx)
    charcont => CastAsCharacterStringType(obj)
    if (associated(charcont)) then
      allocate (character(len=charcont%strlen()) :: string)
      string(:) = charcont
    else
      string = ''
    end if
  end function GetStringFromList

end module StringListModule
