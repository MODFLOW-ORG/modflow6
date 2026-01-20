!> @brief
!! -- Module to manage unit numbers.  Allows for multiple unit numbers
!! -- assigned to a single package type, as shown below.
!! --    row(i)   cunit(i) iunit(i)%nval   iunit(i)%iunit  iunit(i)%ipos
!! --         1      BCF6              1           (1000)            (1)
!! --         2       WEL              3 (1001,1003,1005)        (2,5,7)
!! --         3       GHB              1           (1002)            (4)
!! --         4       EVT              2      (1004,1006)         (6,10)
!! --         5       RIV              0               ()             ()
!! --       ...
!<

module IunitModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENFTYPE
  use SimModule, only: store_error, store_error_filename
  implicit none
  private
  public :: IunitType

  type :: IunitRowType
    integer(I4B) :: nval = 0
    integer(I4B), allocatable, dimension(:) :: iunit ! unit numbers for this row
    integer(I4B), allocatable, dimension(:) :: ipos ! position in the input files character array
  end type IunitRowType

  type :: IunitType
    integer(I4B) :: niunit = 0
    character(len=LENFTYPE), allocatable, dimension(:) :: cunit
    type(IunitRowType), allocatable, dimension(:) :: iunit
  contains
    procedure :: init
    procedure :: addfile
    procedure :: getunitnumber
  end type IunitType

contains

  !> @brief Allocate the cunit and iunit entries of this object, and copy cunit
  !! into the object
  !<
  subroutine init(this, niunit, cunit)
    ! -- dummy
    class(IunitType), intent(inout) :: this
    integer(I4B), intent(in) :: niunit
    character(len=*), dimension(niunit), intent(in) :: cunit
    ! -- local
    integer(I4B) :: i
    !
    allocate (this%cunit(niunit))
    allocate (this%iunit(niunit))
    this%niunit = niunit
    do i = 1, niunit
      this%cunit(i) = cunit(i)
    end do
  end subroutine init

  !> @brief Add an ftyp and unit number
  !!
  !! Find the row for the ftyp and store another iunit value.
  !<
  subroutine addfile(this, ftyp, iunit, ipos, namefilename)
    ! -- dummy
    class(IunitType), intent(inout) :: this
    character(len=*), intent(in) :: ftyp
    integer(I4B), intent(in) :: iunit
    integer(I4B), intent(in) :: ipos
    character(len=*), intent(in) :: namefilename
    ! -- local
    character(len=LINELENGTH) :: errmsg
    integer(I4B), allocatable, dimension(:) :: itemp
    integer(I4B) :: i, irow
    !
    ! -- Find the row containing ftyp
    irow = 0
    do i = 1, this%niunit
      if (this%cunit(i) == ftyp) then
        irow = i
        exit
      end if
    end do
    if (irow == 0) then
      write (errmsg, '(a,a)') 'Package type not supported: ', ftyp
      call store_error(errmsg)
      call store_error_filename(namefilename, terminate=.TRUE.)
    end if
    !
    ! -- Store the iunit number for this ftyp
    if (this%iunit(irow)%nval == 0) then
      allocate (this%iunit(irow)%iunit(1))
      allocate (this%iunit(irow)%ipos(1))
      this%iunit(irow)%nval = 1
    else
      !
      ! -- Increase size of iunit
      allocate (itemp(this%iunit(irow)%nval))
      itemp(:) = this%iunit(irow)%iunit(:)
      deallocate (this%iunit(irow)%iunit)
      this%iunit(irow)%nval = this%iunit(irow)%nval + 1
      allocate (this%iunit(irow)%iunit(this%iunit(irow)%nval))
      this%iunit(irow)%iunit(1:this%iunit(irow)%nval - 1) = itemp(:)
      !
      ! -- Increase size of ipos
      itemp(:) = this%iunit(irow)%ipos(:)
      deallocate (this%iunit(irow)%ipos)
      allocate (this%iunit(irow)%ipos(this%iunit(irow)%nval))
      this%iunit(irow)%ipos(1:this%iunit(irow)%nval - 1) = itemp(:)
      !
      ! -- Cleanup temp
      deallocate (itemp)
    end if
    this%iunit(irow)%iunit(this%iunit(irow)%nval) = iunit
    this%iunit(irow)%ipos(this%iunit(irow)%nval) = ipos
  end subroutine

  !> @brief Get the last unit number for type ftyp or return 0 for iunit.
  !!
  !! If iremove is 1, then remove this unit number.  Similar to a list.pop().
  !<
  subroutine getunitnumber(this, ftyp, iunit, iremove)
    ! -- dummy
    class(IunitType), intent(inout) :: this
    character(len=*), intent(in) :: ftyp
    integer(I4B), intent(inout) :: iunit
    integer(I4B), intent(in) :: iremove
    ! -- local
    integer(I4B) :: i, irow, nval
    !
    ! -- Find the row
    irow = 0
    do i = 1, this%niunit
      if (this%cunit(i) == ftyp) then
        irow = i
        exit
      end if
    end do
    !
    ! -- Find the unit number
    iunit = 0
    if (irow > 0) then
      nval = this%iunit(irow)%nval
      if (nval > 0) then
        iunit = this%iunit(irow)%iunit(nval)
        if (iremove > 0) then
          this%iunit(irow)%iunit(nval) = 0
          this%iunit(irow)%nval = nval - 1
        end if
      else
        iunit = 0
      end if
    end if
  end subroutine getunitnumber

end module IunitModule
