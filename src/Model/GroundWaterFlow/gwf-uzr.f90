module GwfUzrModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENVARNAME, DHALF, DHNOFLO, LINELENGTH
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value
  use SimModule, only: store_error
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use GwfNpfModule, only: GwfNpfType
  use GwfNpfExtModule, only: GwfNpfExtType
  use GwfStoModule, only: GwfStoType
  use GwfStoExtModule, only: GwfStoExtType
  use UzrFlowModule, only: UzrFlowType
  use UzrStorageModule, only: UzrStorageType
  use GwfUzrInputModule, only: GwfUzrParamFoundType
  use UzrSoilModelModule, only: SoilModelType
  use UzrSoilModelFactoryModule, only: init_soil_models, &
                                       get_soil_model, &
                                       destroy_soil_models, &
                                       soil_model_name, &
                                       BROOKS_COREY, HAVERKAMP, VANGENUCHTEN
  implicit none
  private

  public :: GwfUzrType, uzr_cr

  type, extends(NumericalPackageType) :: GwfUzrType
    integer(I4B), pointer :: ioutsat => null() !< output unit for writing saturation to binary file (0 when inactive)
    integer(I4B), pointer :: ioutphead => null() !< output unit for writing pressure head to binary file (0 when inactive)
    integer(I4B), pointer :: storage_scheme => null() !< 0 = default, 1 = chord slope, 2 = mod. Picard
    integer(I4B), pointer :: soil_model_id => null() !< 0 = default, 1 = Brooks-Corey, 2 = Haverkamp, 3 = Van Genuchten
    class(SoilModelType), pointer :: soil_model => null() !< the soil model
    integer(I4B), pointer :: soil_model_kr_id => null() !< separate soil model id for relative permeability (development)
    class(SoilModelType), pointer :: soil_model_kr => null() !< alternative soil model for relative permeability (development)
    integer(I4B), pointer :: kr_averaging => null() !< 0 = default, 1 = geometric, 2 = arithm. mean, 3 = upstream
    integer(I4B), dimension(:), pointer, contiguous :: iunsat => null() !< 0 = standard node, 1 = unsaturated (Richards) node
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< the volumetric fraction of the pore space
    real(DP), dimension(:), pointer, contiguous :: sat_res => null() !< residual (also irreducible) saturation
    real(DP), dimension(:), pointer, contiguous :: pressure_head => null() !< pressure head
    real(DP), dimension(:), pointer, contiguous :: saturation => null() !< water saturation
    class(UzrFlowType), pointer :: uzr_flow => null() !< the NPF flow extension
    class(UzrStorageType), pointer :: uzr_sto => null() !< the STO storage extension
  contains
    procedure :: uzr_df
    procedure :: uzr_ar
    procedure :: uzr_cq
    procedure :: uzr_ot_dv
    procedure :: uzr_da
    procedure :: allocate_scalars
    ! private
    procedure, private :: source_options
    procedure, private :: source_griddata
    procedure, private :: check_griddata
  end type GwfUzrType

contains

  !> @brief Create the UZR package object
  !<
  subroutine uzr_cr(uzr_obj, name_model, input_mempath, inunit, iout)
    type(GwfUzrType), pointer, intent(inout) :: uzr_obj !< the newly created UZR package
    character(len=*), intent(in) :: name_model !< the model name
    character(len=*), intent(in) :: input_mempath !< the memory path of IDM
    integer(I4B), intent(in) :: inunit !< the input file unit for UZR
    integer(I4B), intent(in) :: iout !< the output file unit for the model
    ! local
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'UZR -- RICHARDS FLOW PACKAGE, VERSION 1, 3/30/2015', &
       &' INPUT READ FROM MEMPATH: ', A, /)"

    allocate (uzr_obj)

    call uzr_obj%set_names(1, name_model, 'UZR', 'UZR', input_mempath)
    call uzr_obj%allocate_scalars()

    uzr_obj%inunit = inunit
    uzr_obj%iout = iout
    uzr_obj%ioutsat = 0
    uzr_obj%ioutphead = 0
    uzr_obj%storage_scheme = -1
    uzr_obj%soil_model_id = -1
    uzr_obj%soil_model_kr_id = -1
    uzr_obj%kr_averaging = -1

    ! print a message when enabled
    if (inunit > 0) then
      write (iout, fmtheader) input_mempath
    end if

  end subroutine uzr_cr

  !> @brief Allocates, sources from input, and initializes uzr
  !<
  subroutine uzr_df(this, dis, npf, sto)
    class(GwfUzrType), intent(inout) :: this !< this instance
    class(DisBaseType), pointer, intent(inout) :: dis !< the model discretization
    type(GwfNpfType), pointer, intent(inout) :: npf !< the NPF package to interact with
    type(GwfStoType), pointer, intent(inout) :: sto !< the STO package to interact with
    ! local
    class(GwfNpfExtType), pointer :: npf_ext
    class(GwfStoExtType), pointer :: sto_ext
    integer(I4B) :: i
    character(len=LINELENGTH) :: errmsg

    this%dis => dis

    ! allocate arrays
    call mem_allocate(this%iunsat, dis%nodes, "IUNSAT", this%memoryPath)
    call mem_allocate(this%porosity, dis%nodes, 'POROSITY', this%memoryPath)
    call mem_allocate(this%sat_res, dis%nodes, 'SATRES', this%memoryPath)
    call mem_allocate(this%pressure_head, dis%nodes, 'PHEAD', this%memoryPath)
    call mem_allocate(this%saturation, dis%nodes, 'SATURATION', this%memoryPath)

    ! init
    do i = 1, dis%nodes
      this%pressure_head(i) = DHNOFLO
      this%saturation(i) = DHNOFLO
    end do

    ! load from idm
    call this%source_options()
    call this%source_griddata()

    ! validate
    do i = 1, dis%nodes
      if (this%iunsat(i) == 0) then
        write (errmsg, '(a,a)') &
          'Activating Richards flow in part of the model domain is not yet&
          & supported. Modify IUNSAT in UZR package for model ', &
          trim(this%name_model)
        call store_error(errmsg, .true.)
      end if
    end do

    ! set up the soil model(s)
    call init_soil_models(this%porosity, this%sat_res, this%dis, &
                          this%memoryPath, this%input_mempath)
    this%soil_model => get_soil_model(this%soil_model_id)
    this%soil_model_kr => this%soil_model
    if (this%soil_model_kr_id /= this%soil_model_id) then
      ! dev option: separate model for kr
      this%soil_model_kr => get_soil_model(this%soil_model_kr_id)
    end if

    ! inject Richards into NPF:
    allocate (this%uzr_flow)
    call this%uzr_flow%initialize(this%iunsat, this%kr_averaging, &
                                  this%soil_model_kr, dis, npf)

    npf_ext => this%uzr_flow
    call npf%set_flow_extension(npf_ext)

    ! inject storage calculation
    allocate (this%uzr_sto)
    call this%uzr_sto%initialize(this%iunsat, this%storage_scheme, &
                                 this%soil_model, dis, sto)

    sto_ext => this%uzr_sto
    call sto%set_storage_extension(sto_ext)

  end subroutine uzr_df

  subroutine allocate_scalars(this)
    class(GwfUzrType) :: this

    call this%NumericalPackageType%allocate_scalars()

    call mem_allocate(this%ioutsat, 'IOUTSAT', this%memoryPath)
    call mem_allocate(this%ioutphead, 'IOUTPHEAD', this%memoryPath)
    call mem_allocate(this%soil_model_id, "SOIL_MODEL", this%memoryPath)
    call mem_allocate(this%soil_model_kr_id, "SOIL_MODEL_KR", this%memoryPath)
    call mem_allocate(this%storage_scheme, "STORAGE_SCHEME", this%memoryPath)
    call mem_allocate(this%kr_averaging, "KR_AVERAGING", this%memoryPath)

  end subroutine allocate_scalars

  subroutine source_options(this)
    use UzrFlowModule, only: kr_averaging_name
    use FeatureFlagsModule, only: developmode
    use InputOutputModule, only: getunit, openfile
    use OpenSpecModule, only: access, form
    class(GwfUzrType), intent(inout) :: this
    ! local
    type(GwfUzrParamFoundType) :: found
    character(len=LENVARNAME), dimension(2) :: scheme_name = &
      &[character(len=LENVARNAME) :: 'CHORD-SLOPE', 'MODIFIED-PICARD']
    character(len=LINELENGTH) :: sat_fname, phead_fname

    write (this%iout, '(1x,a)') 'Setting UZR options'

    this%soil_model_id = 0
    call mem_set_value(this%soil_model_id, 'SOIL_MODEL', this%input_mempath, &
                       soil_model_name, found%soil_model)
    if (this%soil_model_id > 0) then
      write (this%iout, '(4x,2a)') 'Soil model set to ', &
        trim(soil_model_name(this%soil_model_id))
    else
      this%soil_model_id = BROOKS_COREY ! default
    end if
    this%soil_model_kr_id = 0
    call mem_set_value(this%soil_model_kr_id, 'MODEL_KR', this%input_mempath, &
                       soil_model_name, found%model_kr)
    if (this%soil_model_kr_id > 0) then
      call developmode('MODEL_KR is a development feature, install the &
        &nightly build or compile from source with IDEVELOPMODE = 1.')
      write (this%iout, '(4x,2a)') &
        'Soil model for relative permeability set to ', &
        trim(soil_model_name(this%soil_model_kr_id))
    else
      this%soil_model_kr_id = this%soil_model_id
    end if

    this%storage_scheme = 0
    if (this%inewton == 0) then
      ! chord-slope and mod. picard only make sense when not newton
      call mem_set_value(this%storage_scheme, 'STORAGE_SCHEME', &
                         this%input_mempath, scheme_name, found%storage_scheme)
    end if
    if (this%storage_scheme > 0) then
      write (this%iout, '(4x,2a)') 'Storage scheme set to ', &
        trim(scheme_name(this%storage_scheme))
    end if

    this%kr_averaging = 0
    call mem_set_value(this%kr_averaging, 'KR_AVERAGING', this%input_mempath, &
                       kr_averaging_name, found%kr_averaging)
    if (this%kr_averaging > 0) then
      write (this%iout, '(4x,2a)') 'Kr weighting set to ', &
        trim(kr_averaging_name(this%kr_averaging))
    end if

    ! set fileout for saturation
    call mem_set_value(sat_fname, 'SATURATIONFILE', this%input_mempath, &
                       found%saturationfile)                       
    if (found%saturationfile) then
      write (this%iout, '(4x,3a)') 'Saturation will be written to ', &
        trim(sat_fname), ' based on HEAD save record in Output Control'
      this%ioutsat = getunit()
      call openfile(this%ioutsat, this%iout, trim(sat_fname), &
                    'DATA(BINARY)', form, access, 'REPLACE')
    end if

    ! set fileout for pressure head
    call mem_set_value(phead_fname, 'PHEADFILE', this%input_mempath, &
                       found%pheadfile)
    
    if (found%pheadfile) then
      write (this%iout, '(4x,3a)') 'Pressure head will be written to ', &
        trim(phead_fname), ' based on HEAD save record in Output Control'
      this%ioutphead = getunit()
      call openfile(this%ioutphead, this%iout, trim(phead_fname), &
                    'DATA(BINARY)', form, access, 'REPLACE')
    end if

    write (this%iout, '(1x,a)') 'End setting UZR options'

  end subroutine source_options

  subroutine source_griddata(this)
    use MemoryManagerExtModule, only: mem_set_value
    class(GwfUzrType), intent(inout) :: this
    ! local
    type(GwfUzrParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map

    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser

    call mem_set_value(this%iunsat, 'IUNSAT', this%input_mempath, &
                       map, found%iunsat)
    call mem_set_value(this%porosity, 'POROSITY', this%input_mempath, &
                       map, found%porosity)
    call mem_set_value(this%sat_res, 'SATRES', this%input_mempath, &
                       map, found%satres)
    call this%check_griddata(found)

  end subroutine source_griddata

  subroutine check_griddata(this, found)
    class(GwfUzrType), intent(inout) :: this
    type(GwfUzrParamFoundType) :: found

  end subroutine check_griddata

  subroutine uzr_ar(this)
    class(GwfUzrType), intent(inout) :: this

  end subroutine uzr_ar

  !< @brief Calculate derived variables
  !>
  subroutine uzr_cq(this, hnew)
    class(GwfUzrType), intent(inout) :: this ! this instance
    real(DP), intent(inout), dimension(:) :: hnew !< the updated head
    ! local
    integer(I4B) :: n
    real(DP) :: z, psi

    do n = 1, this%dis%nodes
      if (this%iunsat(n) == 0) cycle

      ! calculate psi
      z = DHALF * (this%dis%bot(n) + this%dis%top(n))
      psi = hnew(n) - z
      ! calculate saturation
      this%pressure_head(n) = psi
      this%saturation(n) = this%soil_model%saturation(psi, n)
    end do

  end subroutine uzr_cq

  !> @brief Save data to binary file
  !<
  subroutine uzr_ot_dv(this, idvfl)
    class(GwfUzrType) :: this !, this instance
    integer(I4B), intent(in) :: idvfl !< when > 0, dep. var. is written
    ! local
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    integer(I4B) :: ibinun
    integer(I4B) :: nvaluesp
    integer(I4B) :: nwidthp

    if (idvfl > 0) then
      if (this%ioutsat /= 0) then
        ibinun = this%ioutsat
        call this%dis%record_array(this%saturation, this%iout, 0, ibinun, &
                                   '      SATURATION', cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, DHNOFLO)
      end if
      if (this%ioutphead /= 0) then
        ibinun = this%ioutphead
        call this%dis%record_array(this%pressure_head, this%iout, 0, ibinun, &
                                   '       PRESSHEAD', cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, DHNOFLO)
      end if
    end if

  end subroutine uzr_ot_dv

  !< @brief Clean up
  !<
  subroutine uzr_da(this)
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    class(GwfUzrType), intent(inout) :: this !< this instance

    call mem_deallocate(this%ioutsat)
    call mem_deallocate(this%ioutphead)
    call mem_deallocate(this%soil_model_id)
    call mem_deallocate(this%soil_model_kr_id)
    call mem_deallocate(this%storage_scheme)
    call mem_deallocate(this%kr_averaging)

    if (this%inunit > 0) then
      call destroy_soil_models()

      call this%uzr_flow%destroy()
      deallocate (this%uzr_flow)

      call this%uzr_sto%destroy()
      deallocate (this%uzr_sto)

      call mem_deallocate(this%iunsat)
      call mem_deallocate(this%porosity)
      call mem_deallocate(this%sat_res)
      call mem_deallocate(this%pressure_head)
      call mem_deallocate(this%saturation)
    end if

    call memorystore_remove(this%name_model, 'UZR', idm_context)

    call this%NumericalPackageType%da()

  end subroutine uzr_da

end module GwfUzrModule
