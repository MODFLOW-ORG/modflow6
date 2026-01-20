module GwfOcModule

  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMODELNAME
  use OutputControlModule, only: OutputControlType
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr

  implicit none
  private
  public GwfOcType, oc_cr

  !> @ brief Output control for GWF
  !!
  !!  Concrete implementation of OutputControlType for the
  !!  GWF Model
  !<
  type, extends(OutputControlType) :: GwfOcType
  contains
    procedure :: oc_ar
  end type GwfOcType

contains

  !> @ brief Create GwfOcType
  !!
  !!  Create by allocating a new GwfOcType object and initializing
  !!  member variables.
  !!
  !<
  subroutine oc_cr(ocobj, name_model, input_mempath, inunit, iout)
    ! -- dummy
    type(GwfOcType), pointer :: ocobj !< GwfOcType object
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

  !> @ brief Allocate and read GwfOcType
  !!
  !!  Setup head and budget as output control variables.
  !!
  !<
  subroutine oc_ar(this, head, dis, dnodata)
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerExtModule, only: mem_set_value
    use GwfOcInputModule, only: GwfOcParamFoundType
    ! -- dummy
    class(GwfOcType) :: this !< GwfOcType object
    real(DP), dimension(:), pointer, contiguous, intent(in) :: head !< model head
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization package
    real(DP), intent(in) :: dnodata !< no data value
    ! -- local
    integer(I4B) :: i, nocdobj, inodata
    type(OutputControlDataType), pointer :: ocdobjptr
    real(DP), dimension(:), pointer, contiguous :: nullvec => null()
    character(len=LINELENGTH) :: headfile
    type(GwfOcParamFoundType) :: found
    !
    ! -- Initialize variables
    inodata = 0
    nocdobj = 2
    allocate (this%ocds(nocdobj))
    do i = 1, nocdobj
      call ocd_cr(ocdobjptr)
      select case (i)
      case (1)
        call ocdobjptr%init_dbl('BUDGET', nullvec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      case (2)
        call ocdobjptr%init_dbl('HEAD', head, dis, 'PRINT LAST ', &
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
      call mem_set_value(headfile, 'HEADFILE', this%input_mempath, &
                         found%headfile)
      if (found%headfile) then
        call this%set_ocfile('HEAD', headfile, this%iout)
      end if
      write (this%iout, '(1x,a)') 'END OF OC OPTIONS'
    end if
  end subroutine oc_ar

end module GwfOcModule
