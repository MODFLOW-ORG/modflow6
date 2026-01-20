!> @brief Channel Flow (OLF) Module
!<
module OlfModule

  use KindModule, only: I4B
  use ConstantsModule, only: LENPACKAGETYPE, LENMEMPATH, LINELENGTH
  use SimModule, only: store_error
  use BaseModelModule, only: BaseModelType
  use ListsModule, only: basemodellist
  use BaseModelModule, only: AddBaseModelToList
  use SwfModule, only: SwfModelType
  use BudgetModule, only: budget_cr

  implicit none

  private
  public :: olf_cr
  public :: OlfModelType
  public :: OLF_NBASEPKG, OLF_NMULTIPKG
  public :: OLF_BASEPKG, OLF_MULTIPKG

  type, extends(SwfModelType) :: OlfModelType
  contains
    procedure :: set_namfile_options
    procedure :: log_namfile_options
  end type OlfModelType

  !> @brief OLF base package array descriptors
  !!
  !! OLF model base package types.  Only listed packages are candidates
  !< for input and these will be loaded in the order specified.
  integer(I4B), parameter :: OLF_NBASEPKG = 7
  character(len=LENPACKAGETYPE), dimension(OLF_NBASEPKG) :: &
    OLF_BASEPKG = ['DIS2D6 ', 'DISV2D6', 'DFW6   ', &
                   'OC6    ', 'IC6    ', 'OBS6   ', &
                   'STO6   ']

  !> @brief OLF multi package array descriptors
  !!
  !! OLF model multi-instance package types.  Only listed packages are
  !< candidates for input and these will be loaded in the order specified.
  integer(I4B), parameter :: OLF_NMULTIPKG = 50
  character(len=LENPACKAGETYPE), dimension(OLF_NMULTIPKG) :: OLF_MULTIPKG
  data OLF_MULTIPKG/'FLW6 ', 'CHD6 ', 'CDB6 ', 'ZDG6 ', 'PCP6 ', & !  5
                    'EVP6 ', '     ', '     ', '     ', '     ', & !  10
                    &40*'     '/ ! 50

  ! size of supported model package arrays
  integer(I4B), parameter :: NIUNIT_OLF = OLF_NBASEPKG + OLF_NMULTIPKG

contains

  !> @brief Create a new overland flow model object
  !<
  subroutine olf_cr(filename, id, modelname)
    ! modules
    ! dummy
    character(len=*), intent(in) :: filename !< input file
    integer(I4B), intent(in) :: id !< consecutive model number listed in mfsim.nam
    character(len=*), intent(in) :: modelname !< name of the model
    ! local
    class(OlfModelType), pointer :: this
    class(BaseModelType), pointer :: model

    ! Allocate a new model
    allocate (this)
    model => this
    call AddBaseModelToList(basemodellist, model)

    ! call parent initialize routine
    call this%initialize('OLF', filename, id, modelname)

    ! set and log namefile options
    call this%set_namfile_options()

    ! Create utility objects
    call budget_cr(this%budget, this%name)

    ! create model packages
    call this%create_packages()

  end subroutine olf_cr

  !> @brief Handle namefile options
  !!
  !! Set pointers to IDM namefile options, then
  !! create the list file and log options.
  !<
  subroutine set_namfile_options(this)
    use SimVariablesModule, only: idm_context
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use OlfNamInputModule, only: OlfNamParamFoundType
    class(OlfModelType) :: this
    type(OlfNamParamFoundType) :: found
    character(len=LENMEMPATH) :: input_mempath
    character(len=LINELENGTH) :: lst_fname

    ! set input model namfile memory path
    input_mempath = create_mem_path(this%name, 'NAM', idm_context)

    ! copy option params from input context
    call mem_set_value(lst_fname, 'LIST', input_mempath, found%list)
    call mem_set_value(this%inewton, 'NEWTON', input_mempath, found%newton)
    call mem_set_value(this%inewtonur, 'UNDER_RELAXATION', input_mempath, &
                       found%under_relaxation)
    call mem_set_value(this%iprpak, 'PRINT_INPUT', input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', input_mempath, &
                       found%print_flows)
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', input_mempath, found%save_flows)

    ! create the list file
    call this%create_lstfile(lst_fname, this%filename, found%list, &
                             'CHANNEL FLOW MODEL (OLF)')

    ! activate save_flows if found
    if (found%save_flows) then
      this%ipakcb = -1
    end if

    ! log set options
    if (this%iout > 0) then
      call this%log_namfile_options(found)
    end if

  end subroutine set_namfile_options

  !> @brief Write model namfile options to list file
  !<
  subroutine log_namfile_options(this, found)
    use OlfNamInputModule, only: OlfNamParamFoundType
    class(OlfModelType) :: this
    type(OlfNamParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'BEGIN NAMEFILE OPTIONS'

    if (found%print_input) then
      write (this%iout, '(4x,a)') 'STRESS PACKAGE INPUT WILL BE PRINTED '// &
        'FOR ALL MODEL STRESS PACKAGES'
    end if

    if (found%print_flows) then
      write (this%iout, '(4x,a)') 'PACKAGE FLOWS WILL BE PRINTED '// &
        'FOR ALL MODEL PACKAGES'
    end if

    if (found%save_flows) then
      write (this%iout, '(4x,a)') &
        'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') 'END NAMEFILE OPTIONS'

  end subroutine log_namfile_options

end module OlfModule
