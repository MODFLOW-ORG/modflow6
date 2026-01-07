!> @brief Model output control.
module OutputControlModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMODELNAME, LENMEMPATH, LINELENGTH
  use SimVariablesModule, only: errmsg
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr
  use InputOutputModule, only: GetUnit, openfile

  implicit none
  private
  public OutputControlType, oc_cr

  !> @ brief Controls model output. Overridden for each model type.
  type OutputControlType
    character(len=LENMEMPATH) :: memoryPath !< path to data stored in the memory manager
    character(len=LENMEMPATH) :: input_mempath !< input memory path
    character(len=LENMODELNAME), pointer :: name_model => null() !< name of the model
    character(len=LINELENGTH), pointer :: input_fname => null() !< input file name
    integer(I4B), pointer :: inunit => null() !< unit number for input file
    integer(I4B), pointer :: iout => null() !< unit number for output file
    integer(I4B), pointer :: ibudcsv => null() !< unit number for budget csv output file
    integer(I4B), pointer :: iperoc => null() !< stress period number for next output control
    integer(I4B), pointer :: iocrep => null() !< output control repeat flag (period 0 step 0)
    type(OutputControlDataType), pointer, contiguous :: ocds(:) => null() !< output control objects
  contains
    procedure :: oc_df
    procedure :: oc_rp
    procedure :: oc_ot
    procedure :: oc_da
    procedure :: allocate_scalars => allocate
    procedure :: source_options
    procedure :: set_ocfile
    procedure :: oc_save
    procedure :: oc_print
    procedure :: oc_save_unit
    procedure :: set_print_flag
  end type OutputControlType

contains

  !> @brief Create a new output control object.
  subroutine oc_cr(oc, name_model, input_mempath, inunit, iout)
    type(OutputControlType), pointer :: oc !< OutputControlType object
    character(len=*), intent(in) :: name_model !< name of the model
    character(len=*), intent(in) :: input_mempath !< input mempath of the package
    integer(I4B), intent(in) :: inunit !< unit number for input
    integer(I4B), intent(in) :: iout !< unit number for output

    allocate (oc)
    call oc%allocate_scalars(name_model, input_mempath)
    oc%inunit = inunit
    oc%iout = iout
  end subroutine oc_cr

  !> @ brief Define the output control type. Placeholder routine.
  subroutine oc_df(this)
    class(OutputControlType) :: this !< this instance
  end subroutine oc_df

  !> @ brief Read period block options and prepare the output control type.
  subroutine oc_rp(this)
    ! modules
    use TdisModule, only: kper
    use SimModule, only: store_error, store_error_filename
    use MemoryManagerModule, only: mem_setptr
    use CharacterStringModule, only: CharacterStringType
    ! dummy
    class(OutputControlType) :: this !< this instance
    ! local
    class(OutputControlDataType), pointer :: ocdobjptr
    integer(I4B), pointer :: iper
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: ocactions
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: rtypes
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: ocsettings
    integer(I4B), pointer :: nlist
    integer(I4B) :: n, ipos
    character(len=LINELENGTH) :: ocaction, rtype, ocsetting
    logical(LGP) :: found_rtype
    ! formats
    character(len=*), parameter :: fmtboc = &
      &"(1X,/1X,'BEGIN READING OUTPUT CONTROL FOR STRESS PERIOD ',I0)"
    character(len=*), parameter :: fmteoc = &
      &"(/,1X,'END READING OUTPUT CONTROL FOR STRESS PERIOD ',I0)"
    character(len=*), parameter :: fmtroc = &
      "(1X,/1X,'OUTPUT CONTROL FOR STRESS PERIOD ',I0, &
      &' IS REPEATED USING SETTINGS FROM A PREVIOUS STRESS PERIOD.')"

    if (this%input_mempath == '') return
    call mem_setptr(iper, 'IPER', this%input_mempath)
    if (iper /= kper) then
      ! previous output control settings are still active
      write (this%iout, fmtroc) kper
      return
    else
      this%iperoc = iper
      write (this%iout, fmtboc) this%iperoc
    end if

    ! Clear io flags
    do ipos = 1, size(this%ocds)
      ocdobjptr => this%ocds(ipos)
      call ocdobjptr%psm%init()
    end do

    call mem_setptr(nlist, 'NBOUND', this%input_mempath)
    call mem_setptr(ocactions, 'OCACTION', this%input_mempath)
    call mem_setptr(rtypes, 'RTYPE', this%input_mempath)
    call mem_setptr(ocsettings, 'SETTING', this%input_mempath)

    do n = 1, nlist
      ocaction = ocactions(n)
      rtype = rtypes(n)
      ocsetting = ocsettings(n)

      found_rtype = .false.
      do ipos = 1, size(this%ocds)
        ocdobjptr => this%ocds(ipos)
        if (rtype == trim(ocdobjptr%cname)) then
          found_rtype = .true.
          call ocdobjptr%psm%set(ocaction, ocsetting, this%iout)
          call ocdobjptr%ocd_rp_check(this%input_fname)
        end if
      end do
      if (.not. found_rtype) then
        call store_error('Input OC period block rtype not found:  "'// &
                         trim(rtype)//'".')
        call store_error_filename(this%input_fname)
      end if
    end do

    write (this%iout, fmteoc) this%iperoc
  end subroutine oc_rp

  !> @ brief Write output.
  !!
  !! Go through each output control data type and print
  !! and/or save data based on user-specified controls.
  !<
  subroutine oc_ot(this, ipflg)
    ! modules
    use TdisModule, only: kstp, endofperiod
    ! dummy
    class(OutputControlType) :: this !< OutputControlType object
    integer(I4B), intent(inout) :: ipflg !< flag indicating if data was printed
    ! local
    integer(I4B) :: ipos
    type(OutputControlDataType), pointer :: ocdobjptr

    ! Clear printout flag(ipflg).  This flag indicates that an array was
    ! printed to the listing file.
    ipflg = 0

    do ipos = 1, size(this%ocds)
      ocdobjptr => this%ocds(ipos)
      call ocdobjptr%ocd_ot(ipflg, kstp, endofperiod, this%iout)
    end do
  end subroutine oc_ot

  !> @ brief Deallocate method for OutputControlType
  !!
  !!  Deallocate member variables.
  !!
  !<
  subroutine oc_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    ! dummy
    class(OutputControlType) :: this !< OutputControlType object
    ! local
    integer(I4B) :: i

    do i = 1, size(this%ocds)
      call this%ocds(i)%ocd_da()
    end do
    deallocate (this%ocds)

    deallocate (this%name_model)
    deallocate (this%input_fname)
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%iperoc)
    call mem_deallocate(this%iocrep)
  end subroutine oc_da

  !> @ brief Allocate variables for the output control object
  subroutine allocate (this, name_model, input_mempath)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    ! dummy
    class(OutputControlType) :: this !< this instance
    character(len=*), intent(in) :: name_model !< name of model
    character(len=*), intent(in) :: input_mempath !< input mempath of the package
    logical(LGP) :: found

    this%memoryPath = create_mem_path(name_model, 'OC')

    allocate (this%name_model)
    allocate (this%input_fname)
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%iperoc, 'IPEROC', this%memoryPath)
    call mem_allocate(this%iocrep, 'IOCREP', this%memoryPath)

    this%name_model = name_model
    this%input_mempath = input_mempath
    this%input_fname = ''
    this%inunit = 0
    this%iout = 0
    this%ibudcsv = 0
    this%iperoc = 0
    this%iocrep = 0

    if (this%input_mempath /= '') then
      call mem_set_value(this%input_fname, 'INPUT_FNAME', &
                         this%input_mempath, found)
    end if
  end subroutine allocate

  !> @ brief Read the output control options block
  subroutine source_options(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    ! dummy
    class(OutputControlType) :: this !< this instance
    type(OutputControlDataType), pointer :: ocdobjptr
    character(len=LINELENGTH) :: budgetfn, budgetcsv
    character(len=LINELENGTH) :: prnfmt, print_format
    logical(LGP) :: found_budcsv, found_budget
    logical(LGP), dimension(4) :: found_format
    integer(I4B), pointer :: columns, width, ndigits
    integer(I4B) :: ipos

    found_format = .false.
    allocate (columns)
    allocate (width)
    allocate (ndigits)

    call mem_set_value(columns, 'COLUMNS', this%input_mempath, &
                       found_format(1))
    call mem_set_value(width, 'WIDTH', this%input_mempath, &
                       found_format(2))
    call mem_set_value(ndigits, 'DIGITS', this%input_mempath, &
                       found_format(3))
    call mem_set_value(prnfmt, 'FORMAT', this%input_mempath, &
                       found_format(4))
    call mem_set_value(budgetcsv, 'BUDGETCSVFILE', this%input_mempath, &
                       found_budcsv)
    call mem_set_value(budgetfn, 'BUDGETFILE', this%input_mempath, &
                       found_budget)

    if (found_budcsv) then
      this%ibudcsv = GetUnit()
      call openfile(this%ibudcsv, this%iout, budgetcsv, 'CSV', &
                    filstat_opt='REPLACE')
    end if

    if (found_budget) then
      call this%set_ocfile('BUDGET', budgetfn, this%iout)
    end if

    if (found_format(1) .and. &
        found_format(2) .and. &
        found_format(3) .and. &
        found_format(4)) then
      write (print_format, '(a,i0,a,i0,a,i0,a)') 'COLUMNS ', columns, &
        ' WIDTH ', width, ' DIGITS ', ndigits, ' '//trim(prnfmt)//' '
      do ipos = 1, size(this%ocds)
        ocdobjptr => this%ocds(ipos)
        if (ocdobjptr%cname /= 'BUDGET') then
          call ocdobjptr%set_prnfmt(print_format, 0)
        end if
      end do
    end if

    deallocate (columns)
    deallocate (width)
    deallocate (ndigits)
  end subroutine source_options

  subroutine set_ocfile(this, cname, ocfile, iout)
    ! modules
    use SimModule, only: store_error, store_error_filename
    ! dummy
    class(OutputControlType) :: this !< OutputControlDataType object
    character(len=*), intent(in) :: cname !< data object cname
    character(len=*), intent(in) :: ocfile !< OC output filename
    integer(I4B), intent(in) :: iout !< Unit number for output
    type(OutputControlDataType), pointer :: ocdobjptr
    integer(I4B) :: ipos
    logical(LGP) :: found
    found = .false.
    do ipos = 1, size(this%ocds)
      ocdobjptr => this%ocds(ipos)
      if (cname == trim(ocdobjptr%cname)) then
        found = .true.
        call ocdobjptr%set_ocfile(ocfile, iout)
      end if
    end do
    if (.not. found) then
      call store_error('OC internal error: oc data type not found for name "'// &
                       trim(cname)//'".')
      call store_error_filename(this%input_fname)
    end if
  end subroutine set_ocfile

  !> @ brief Determine if it is time to save.
  logical function oc_save(this, cname)
    ! modules
    use TdisModule, only: kstp, endofperiod
    ! dummy
    class(OutputControlType) :: this !< OutputControlType object
    character(len=*), intent(in) :: cname !< character string for data name
    ! local
    integer(I4B) :: ipos
    logical(LGP) :: found
    class(OutputControlDataType), pointer :: ocdobjptr
    !
    oc_save = .false.
    found = .false.
    do ipos = 1, size(this%ocds)
      ocdobjptr => this%ocds(ipos)
      if (cname == trim(ocdobjptr%cname)) then
        found = .true.
        exit
      end if
    end do
    if (found) then
      oc_save = ocdobjptr%psm%should_save(kstp, endofperiod)
    end if
  end function oc_save

  !> @ brief Determine if it is time to print.
  logical function oc_print(this, cname)
    ! modules
    use TdisModule, only: kstp, endofperiod
    ! dummy
    class(OutputControlType) :: this !< OutputControlType object
    character(len=*), intent(in) :: cname !< character string for data name
    ! local
    integer(I4B) :: ipos
    logical(LGP) :: found
    class(OutputControlDataType), pointer :: ocdobjptr

    oc_print = .false.
    found = .false.
    do ipos = 1, size(this%ocds)
      ocdobjptr => this%ocds(ipos)
      if (cname == trim(ocdobjptr%cname)) then
        found = .true.
        exit
      end if
    end do
    if (found) then
      oc_print = ocdobjptr%psm%should_print(kstp, endofperiod)
    end if
  end function oc_print

  !> @ brief Determine unit number for saving
  function oc_save_unit(this, cname)
    ! -- modules
    ! -- return
    integer(I4B) :: oc_save_unit
    ! -- dummy
    class(OutputControlType) :: this !< OutputControlType object
    character(len=*), intent(in) :: cname !< character string for data name
    ! -- local
    integer(I4B) :: ipos
    logical(LGP) :: found
    class(OutputControlDataType), pointer :: ocdobjptr
    !
    oc_save_unit = 0
    found = .false.
    do ipos = 1, size(this%ocds)
      ocdobjptr => this%ocds(ipos)
      if (cname == trim(ocdobjptr%cname)) then
        found = .true.
        exit
      end if
    end do
    if (found) then
      oc_save_unit = ocdobjptr%idataun
    end if
  end function oc_save_unit

  !> @ brief Set the print flag based on convergence and other parameters
  function set_print_flag(this, cname, icnvg, endofperiod) result(iprint_flag)
    ! -- modules
    use SimVariablesModule, only: isimcontinue
    ! -- return
    integer(I4B) :: iprint_flag
    ! -- dummy
    class(OutputControlType) :: this !< OutputControlType object
    character(len=*), intent(in) :: cname !< character string for data name
    integer(I4B), intent(in) :: icnvg !< convergence flag
    logical, intent(in) :: endofperiod !< end of period logical flag
    ! -- local
    !
    ! -- default is to not print
    iprint_flag = 0
    !
    ! -- if the output control file indicates that cname should be printed
    if (this%oc_print(cname)) iprint_flag = 1
    !
    ! -- if it is not a CONTINUE run, then set to print if not converged
    if (isimcontinue == 0) then
      if (icnvg == 0) iprint_flag = 1
    end if
    !
    ! -- if it's the end of the period, then set flag to print
    if (endofperiod) iprint_flag = 1
  end function set_print_flag

end module OutputControlModule
