module SwfOcModule

  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B
  use OutputControlModule, only: OutputControlType
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr

  implicit none
  private
  public SwfOcType, oc_cr

  !> @ brief Output control
  !!
  !!  Concrete implementation of OutputControlType
  !<
  type, extends(OutputControlType) :: SwfOcType
  contains
    procedure :: oc_ar
  end type SwfOcType

contains

  !> @ brief Create SwfOcType
  !!
  !!  Create by allocating a new SwfOcType object and initializing
  !!  member variables.
  !!
  !<
  subroutine oc_cr(ocobj, name_model, input_mempath, inunit, iout)
    ! -- dummy
    type(SwfOcType), pointer :: ocobj !< SwfOcType object
    character(len=*), intent(in) :: name_model !< name of the model
    character(len=*), intent(in) :: input_mempath !< input mempath of the package
    integer(I4B), intent(in) :: inunit !< unit number for input
    integer(I4B), intent(in) :: iout !< unit number for output
    !
    ! -- Create the object
    allocate (ocobj)
    !
    ! -- Allocate scalars
    call ocobj%allocate_scalars(name_model, input_mempath)
    !
    ! -- Save unit numbers
    ocobj%inunit = inunit
    ocobj%iout = iout
  end subroutine oc_cr

  !> @ brief Allocate and read SwfOcType
  !!
  !!  Setup head and budget as output control variables.
  !!
  !<
  subroutine oc_ar(this, datavec, dis, dnodata)
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerExtModule, only: mem_set_value
    use SwfOcInputModule, only: SwfOcParamFoundType
    ! -- dummy
    class(SwfOcType) :: this !< SwfOcType object
    real(DP), dimension(:), pointer, contiguous, intent(in) :: datavec !< data vector
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization package
    real(DP), intent(in) :: dnodata !< no data value
    ! -- local
    integer(I4B) :: i, nocdobj, inodata
    type(OutputControlDataType), pointer :: ocdobjptr
    real(DP), dimension(:), pointer, contiguous :: nullvec => null()
    character(len=LINELENGTH) :: stagefile, qoutflowfile
    type(SwfOcParamFoundType) :: found
    !
    ! -- Initialize variables
    inodata = 0
    nocdobj = 3
    allocate (this%ocds(nocdobj))
    do i = 1, nocdobj
      call ocd_cr(ocdobjptr)
      select case (i)
      case (1)
        call ocdobjptr%init_dbl('BUDGET', nullvec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      case (2)
        call ocdobjptr%init_dbl('STAGE', datavec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      case (3)
        call ocdobjptr%init_dbl('QOUTFLOW', datavec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      end select
      this%ocds(i) = ocdobjptr
      deallocate (ocdobjptr)
    end do
    !
    ! -- Read options or set defaults if this package not on
    if (this%input_mempath /= '') then
      write (this%iout, '(/,1x,a,/)') 'PROCESSING OC OPTIONS'
      call this%source_options()
      call mem_set_value(qoutflowfile, 'QOUTFLOWFILE', this%input_mempath, &
                         found%qoutflowfile)
      call mem_set_value(stagefile, 'STAGEFILE', this%input_mempath, &
                         found%stagefile)
      if (found%qoutflowfile) then
        call this%set_ocfile('QOUTFLOW', qoutflowfile, this%iout)
      end if
      if (found%stagefile) then
        call this%set_ocfile('STAGE', stagefile, this%iout)
      end if
      write (this%iout, '(1x,a)') 'END OF OC OPTIONS'
    end if
  end subroutine oc_ar

end module SwfOcModule
