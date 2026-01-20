!> @brief This module contains the MvrModule Module
!!
!! This module contains the code for the low-level MvrType
!! object.
!!
!<
module MvrModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMODELNAME, LINELENGTH, LENBUDTXT, &
                             LENAUXNAME, LENBOUNDNAME, DZERO, DONE, &
                             LENMEMPATH
  use SimVariablesModule, only: errmsg
  use PackageMoverModule, only: PackageMoverType

  implicit none
  private
  public :: MvrType

  character(len=12), dimension(4) :: mvrtypes = &
    &[character(len=12) :: 'FACTOR', 'EXCESS', 'THRESHOLD', 'UPTO']

  !> @brief Derived type for MvrType
  !!
  !! This derived type contains information and methods for
  !! moving water between packages.
  !!
  !<
  type MvrType
    character(len=LENMEMPATH) :: mem_path_src = '' !< provider package name
    character(len=LENMEMPATH) :: mem_path_tgt = '' !< receiver package name
    integer(I4B), pointer :: iRchNrSrc => null() !< provider reach number
    integer(I4B) :: iRchNrSrcMapped !< mapped provider reach number (currently for lake outlet)
    integer(I4B), pointer :: iRchNrTgt => null() !< receiver reach number
    integer(I4B), pointer :: imvrtype => null() !< mover type (1, 2, 3, 4) corresponds to mvrtypes
    real(DP), pointer :: value => null() !< factor or rate depending on mvrtype
    logical(LGP) :: is_provider_active = .true.
    logical(LGP) :: is_receiver_active = .true.
    real(DP) :: qpactual = DZERO !< rate provided to the receiver
    real(DP) :: qavailable = DZERO !< rate available at time of providing
    real(DP), pointer :: qtformvr_ptr => null() !< pointer to total available flow (qtformvr)
    real(DP), pointer :: qformvr_ptr => null() !< pointer to available flow after consumed (qformvr)
    real(DP), pointer :: qtomvr_ptr => null() !< pointer to provider flow rate (qtomvr)
    real(DP), pointer :: qfrommvr_ptr => null() !< pointer to receiver flow rate (qfrommvr)
  contains
    procedure :: set_values
    procedure :: prepare
    procedure :: echo
    procedure :: advance
    procedure :: update_provider
    procedure :: update_receiver
    procedure :: qrcalc
    procedure :: writeflow
  end type MvrType

contains

  !> @ brief Set values from input data
  !!
  !! Set values and pointers for mover object.
  !!
  !<
  subroutine set_values(this, mname1, pname1, id1, mname2, pname2, &
                        id2, imvrtype, value)
    use MemoryHelperModule, only: create_mem_path
    class(MvrType) :: this
    character(len=*), intent(in) :: mname1
    character(len=*), intent(in) :: pname1
    integer(I4B), intent(in), target :: id1
    character(len=*), intent(in) :: mname2
    character(len=*), intent(in) :: pname2
    integer(I4B), intent(in), target :: id2
    integer(I4B), intent(in), target :: imvrtype
    real(DP), intent(in), target :: value

    this%mem_path_src = create_mem_path(mname1, pname1)
    this%iRchNrSrc => id1
    this%mem_path_tgt = create_mem_path(mname2, pname2)
    this%iRchNrTgt => id2
    this%imvrtype => imvrtype
    this%value => value

    ! to be set later
    this%iRchNrSrcMapped = -1
  end subroutine set_values

  !> @ brief Prepare object
  !!
  !! Set values and pointers for mover object.
  !! pckMemPaths is an array of strings which are the memory paths for those
  !! packages. They are composed of model names and package names. The mover
  !! entries must be in pckMemPaths, or this routine will terminate with an error.
  !<
  subroutine prepare(this, inunit, pckMemPaths, pakmovers)
    ! -- modules
    use SimModule, only: store_error, store_error_unit, count_errors
    ! -- dummy
    class(MvrType) :: this !< MvrType object
    integer(I4B), intent(in) :: inunit !< input file unit number
    character(len=LENMEMPATH), &
      dimension(:), pointer, contiguous :: pckMemPaths !< array of strings
    type(PackageMoverType), dimension(:), pointer, contiguous :: pakmovers !< Array of package mover objects
    ! -- local
    real(DP), dimension(:), pointer, contiguous :: temp_ptr => null()
    logical :: found
    integer(I4B) :: i
    integer(I4B) :: ipakloc1, ipakloc2
    !
    ! -- Check to make sure provider and receiver are not the same
    if (this%mem_path_src == this%mem_path_tgt .and. &
        this%iRchNrSrc == this%iRchNrTgt) then
      call store_error('Provider and receiver are the same: '// &
                       trim(this%mem_path_src)//' : '//trim(this%mem_path_tgt))
      call store_error_unit(inunit)
    end if
    !
    ! -- Check to make sure pname1 and pname2 are both listed in pckMemPaths
    !    pname1 is the provider package; pname2 is the receiver package
    found = .false.
    ipakloc1 = 0
    do i = 1, size(pckMemPaths)
      if (this%mem_path_src == pckMemPaths(i)) then
        found = .true.
        ipakloc1 = i
        exit
      end if
    end do
    if (.not. found) then
      call store_error('Mover capability not activated in '//this%mem_path_src)
      call store_error('Add "MOVER" keyword to package options block.')
    end if
    found = .false.
    ipakloc2 = 0
    do i = 1, size(pckMemPaths)
      if (this%mem_path_tgt == pckMemPaths(i)) then
        found = .true.
        ipakloc2 = i
        exit
      end if
    end do
    if (.not. found) then
      call store_error('Mover capability not activated in '//this%mem_path_tgt)
      call store_error('Add "MOVER" keyword to package options block.')
    end if
    if (count_errors() > 0) then
      call store_error_unit(inunit)
    end if

    if (this%is_provider_active) then
      !
      ! -- Set pointer to QTOMVR array in the provider boundary package
      temp_ptr => pakmovers(ipakloc1)%qtomvr
      if (this%iRchNrSrc < 1 .or. this%iRchNrSrc > size(temp_ptr)) then
        call store_error('Provider ID < 1 or greater than package size ')
        write (errmsg, '(a,i0,a,i0)') 'Provider ID = ', this%iRchNrSrc, &
          '; Package size = ', size(temp_ptr)
        call store_error(trim(errmsg))
        call store_error_unit(inunit)
      end if
      this%qtomvr_ptr => temp_ptr(this%iRchNrSrc)
      !
      ! -- Set pointer to QFORMVR array in the provider boundary package
      temp_ptr => pakmovers(ipakloc1)%qformvr
      this%qformvr_ptr => temp_ptr(this%iRchNrSrc)
      !
      ! -- Set pointer to QTFORMVR array in the provider boundary package
      temp_ptr => pakmovers(ipakloc1)%qtformvr
      this%qtformvr_ptr => temp_ptr(this%iRchNrSrc)
    end if

    if (this%is_receiver_active) then
      !
      ! -- Set pointer to QFROMMVR array in the receiver boundary package
      temp_ptr => pakmovers(ipakloc2)%qfrommvr
      if (this%iRchNrTgt < 1 .or. this%iRchNrTgt > size(temp_ptr)) then
        call store_error('Receiver ID < 1 or greater than package size ')
        write (errmsg, '(a,i0,a,i0)') 'Receiver ID = ', this%iRchNrTgt, &
          '; package size = ', size(temp_ptr)
        call store_error(trim(errmsg))
        call store_error_unit(inunit)
      end if
      this%qfrommvr_ptr => temp_ptr(this%iRchNrTgt)
    end if
  end subroutine prepare

  !> @ brief Echo data to list file
  !!
  !! Write mover values to output file.
  !!
  !<
  subroutine echo(this, iout)
    ! -- modules
    ! -- dummy
    class(MvrType) :: this !< MvrType
    integer(I4B), intent(in) :: iout !< unit number for output file
    ! -- local
    !
    write (iout, '(4x, a, a, a, i0)') 'FROM PACKAGE: ', trim(this%mem_path_src), &
      ' FROM ID: ', this%iRchNrSrc
    write (iout, '(4x, a, a, a, i0)') 'TO PACKAGE: ', trim(this%mem_path_tgt), &
      ' TO ID: ', this%iRchNrTgt
    write (iout, '(4x, a, a, a, 1pg15.6,/)') 'MOVER TYPE: ', &
      trim(mvrtypes(this%imvrtype)), ' ', this%value
  end subroutine echo

  !> @ brief Advance
  !!
  !! Advance mover object.  Does nothing now.
  !!
  !<
  subroutine advance(this)
    ! -- modules
    ! -- dummy
    class(MvrType) :: this
    ! -- local
  end subroutine advance

  !> @ brief Formulate coefficients
  !!
  !! Make mover calculations for provider.
  !!
  !<
  subroutine update_provider(this)
    ! -- modules
    ! -- dummy
    class(MvrType) :: this !< MvrType
    ! -- local
    real(DP) :: qavailable, qtanew, qpactual
    !
    ! -- Set qa and this%qavailable equal to available water in package (qtomvr)
    qavailable = this%qformvr_ptr
    qtanew = this%qtformvr_ptr
    this%qavailable = qavailable
    !
    ! -- Using the mover rules, calculate how much of the available water will
    !    be provided from the mover to the receiver.
    qpactual = this%qrcalc(qavailable, qtanew)
    !
    ! -- Store qpactual
    this%qpactual = qpactual
    !
    ! -- Add the calculated qpactual term directly into the provider package
    !    qtomvr array.
    this%qtomvr_ptr = this%qtomvr_ptr + qpactual
    !
    ! -- Reduce the amount of water that is available in the provider package
    !    qformvr array.
    this%qformvr_ptr = this%qformvr_ptr - qpactual
  end subroutine update_provider

  !> @ brief Formulate coefficients
  !!
  !! Make mover calculations for receiver.
  !!
  !<
  subroutine update_receiver(this)
    class(MvrType) :: this !< MvrType
    ! -- Add the calculated qpactual term directly into the receiver package
    !    qfrommvr array.
    this%qfrommvr_ptr = this%qfrommvr_ptr + this%qpactual
  end subroutine update_receiver

  !> @ brief Flow to receiver
  !!
  !! Calculate the rate of water provided to receiver.
  !!
  !<
  function qrcalc(this, qa, qta) result(qr)
    ! -- modules
    ! -- return
    real(DP) :: qr
    ! -- dummy
    class(MvrType) :: this !< MvrType
    real(DP), intent(in) :: qa !< actual flow
    real(DP), intent(in) :: qta !< total available flow
    ! -- local
    ! -- Using the mover rules, calculate how much of the available water will
    !    go to the receiver.
    qr = DZERO
    ! -- Calculate qr
    select case (this%imvrtype)
    case (1)
      ! -- FACTOR uses total available to make calculation, and then
      !    limits qr by consumed available
      if (qta > DZERO) qr = qta * this%value
      qr = min(qr, qa)
    case (2)
      ! -- EXCESS
      if (qa > this%value) then
        qr = qa - this%value
      else
        qr = DZERO
      end if
    case (3)
      ! -- THRESHOLD
      if (this%value > qa) then
        qr = DZERO
      else
        qr = this%value
      end if
    case (4)
      ! -- UPTO
      if (qa > this%value) then
        qr = this%value
      else
        qr = qa
      end if
    end select
  end function qrcalc

  !> @ brief Write flow
  !!
  !! Write a line of output for this mover object.
  !!
  !<
  subroutine writeflow(this, iout)
    ! -- modules
    ! -- dummy
    class(MvrType) :: this !< MvrType
    integer(I4B), intent(in) :: iout !< output file unit number
    ! -- local
    character(len=*), parameter :: fmt = &
      "(1x, a, ' ID ', i0, ' AVAILABLE ', 1(1pg15.6), &
      &' PROVIDED ', 1(1pg15.6), ' TO ', a, ' ID ', i0)"
    !
    write (iout, fmt) trim(this%mem_path_src), this%iRchNrSrc, this%qavailable, &
      this%qpactual, trim(this%mem_path_tgt), this%iRchNrTgt
  end subroutine writeflow

end module MvrModule

