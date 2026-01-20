!> @brief This module contains the SourceLoadModule
!!
!! This module contains the routines needed to generate
!! a loader object for an input source and routines
!! that distribute processing to a particular source.
!!
!<
module SourceLoadModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg, iout
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENMODELNAME, LENFTYPE, &
                             LENPACKAGETYPE, LENPACKAGENAME
  use SimModule, only: store_error, store_error_filename
  use ModflowInputModule, only: ModflowInputType, getModflowInput
  use NCFileVarsModule, only: NCFileVarsType

  implicit none
  private
  public :: create_input_loader
  public :: open_source_file
  public :: load_modelnam, load_simnam, load_simtdis
  public :: remote_model_ndim
  public :: export_cr, export_da
  public :: export_post_prepare, export_post_step
  public :: nc_close
  public :: netcdf_context

contains

  !> @brief factory function to create and setup model package static loader
  !<
  function create_input_loader(component_type, subcomponent_type, &
                               component_name, subcomponent_name, input_type, &
                               load_scope, input_fname, component_fname, &
                               nc_vars) &
    result(loader)
    use SourceCommonModule, only: package_source_type, idm_subcomponent_name
    use InputLoadTypeModule, only: StaticPkgLoadBaseType
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: subcomponent_name
    character(len=*), intent(in) :: input_type
    character(len=*), intent(in) :: load_scope
    character(len=*), intent(in) :: input_fname
    character(len=*), intent(in) :: component_fname
    type(NCFileVarsType), pointer, optional, intent(in) :: nc_vars
    class(StaticPkgLoadBaseType), pointer :: loader
    type(ModflowInputType) :: mf6_input
    character(len=LENPACKAGENAME) :: source_type
    character(len=LENPACKAGENAME) :: sc_name

    ! set subcomponent name
    sc_name = idm_subcomponent_name(component_type, subcomponent_type, &
                                    subcomponent_name)
    ! create description of input
    mf6_input = getModflowInput(input_type, component_type, subcomponent_type, &
                                component_name, sc_name, load_scope, &
                                input_fname)
    ! set package source
    source_type = package_source_type(input_fname)

    ! set source loader for model package
    loader => package_loader(source_type)

    ! initialize loader
    call loader%init(mf6_input, component_name, component_fname, input_fname)

    ! initialize loader netcdf variables data structure
    if (present(nc_vars)) then
      call nc_vars%create_varlists(component_name, sc_name, loader%nc_vars)
    else
      call loader%nc_vars%init(component_name)
    end if
  end function create_input_loader

  !> @brief allocate source model package static loader
  !<
  function package_loader(source_type) result(loader)
    use InputLoadTypeModule, only: StaticPkgLoadBaseType
    use IdmMf6FileModule, only: Mf6FileStaticPkgLoadType
    character(len=*), intent(inout) :: source_type
    class(Mf6FileStaticPkgLoadType), pointer :: mf6file_loader
    class(StaticPkgLoadBaseType), pointer :: loader

    ! initialize
    nullify (loader)

    ! allocate derived object
    select case (source_type)
    case ('MF6FILE')
      allocate (mf6file_loader)
      loader => mf6file_loader
      allocate (loader%nc_vars)
    case default
      write (errmsg, '(a)') &
        'Simulation package input source type "'//trim(source_type)// &
        '" not currently supported.'
      call store_error(errmsg, .true.)
    end select
  end function package_loader

  function open_source_file(pkgtype, filename, modelfname, iout) result(fd)
    use SourceCommonModule, only: package_source_type
    use IdmMf6FileModule, only: open_mf6file
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: modelfname
    integer(I4B), intent(in) :: iout
    integer(I4B) :: fd
    character(len=LENPACKAGENAME) :: source_type

    ! initialize
    fd = 0

    ! set source type
    source_type = package_source_type(filename)
    !
    select case (source_type)
    case ('MF6FILE')
      fd = open_mf6file(pkgtype, filename, modelfname, iout)
    case default
    end select
  end function open_source_file

  subroutine load_modelnam(mtype, mfname, mname, iout)
    use SimVariablesModule, only: simfile
    use SourceCommonModule, only: package_source_type, idm_component_type
    use IdmMf6FileModule, only: input_load
    character(len=*), intent(in) :: mtype
    character(len=*), intent(in) :: mfname
    character(len=*), intent(in) :: mname
    integer(I4B), intent(in) :: iout
    type(ModflowInputType) :: mf6_input
    character(len=LENPACKAGENAME) :: source_type

    ! set source type
    source_type = package_source_type(mfname)

    ! create description of input
    mf6_input = getModflowInput(mtype, idm_component_type(mtype), 'NAM', &
                                mname, 'NAM', 'SIM', mfname)
    select case (source_type)
    case ('MF6FILE')
      call input_load(mfname, mf6_input, simfile, iout)
    case default
    end select
  end subroutine load_modelnam

  subroutine load_simnam()
    use SimVariablesModule, only: simfile, iout
    use MemoryManagerModule, only: mem_setptr
    use MessageModule, only: write_message
    use IdmMf6FileModule, only: input_load
    use SourceCommonModule, only: filein_fname
    type(ModflowInputType) :: mf6_input, hpc_input
    character(len=LINELENGTH) :: hpc6_filename
    character(len=LINELENGTH) :: line
    logical(LGP) :: lexist

    ! load mfsim.nam if it exists
    inquire (file=trim(adjustl(simfile)), exist=lexist)

    if (lexist) then
      ! write name of namfile to stdout
      write (line, '(2(1x,a))') 'Using Simulation name file:', &
        trim(adjustl(simfile))
      call write_message(line, skipafter=1)
      ! create description of input
      mf6_input = getModflowInput('NAM6', 'SIM', 'NAM', 'SIM', 'NAM', 'ROOT', &
                                  simfile)
      ! open namfile and load to input context
      call input_load(simfile, mf6_input, simfile, iout)
      ! load optional HPC configuration file
      if (filein_fname(hpc6_filename, 'HPC6_FILENAME', mf6_input%mempath, &
                       simfile)) then
        hpc_input = getModflowInput('HPC6', 'UTL', 'HPC', 'UTL', 'HPC', 'SIM')
        call input_load(hpc6_filename, hpc_input, simfile, iout)
      end if
    end if
  end subroutine load_simnam

  subroutine load_simtdis()
    use SimVariablesModule, only: simfile, iout
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr
    use SimVariablesModule, only: idm_context
    use SourceCommonModule, only: package_source_type
    use IdmMf6FileModule, only: input_load
    character(len=LENMEMPATH) :: input_mempath
    type(ModflowInputType) :: mf6_input
    character(len=LENPACKAGENAME) :: source_type
    character(len=:), pointer :: tdis6
    logical(LGP) :: lexist

    ! set input memory path
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)

    ! set pointers to input context timing params
    call mem_setptr(tdis6, 'TDIS6', input_mempath)

    ! create timing
    if (tdis6 /= '') then
      ! set source type
      source_type = package_source_type(tdis6)
      select case (source_type)
      case ('MF6FILE')
        inquire (file=trim(adjustl(tdis6)), exist=lexist)
        if (lexist) then
          ! create description of input
          mf6_input = getModflowInput('TDIS6', 'SIM', 'TDIS', &
                                      'SIM', 'TDIS', 'SIM', simfile)
          ! open namfile and load to input context
          call input_load(tdis6, mf6_input, simfile, iout)
        else
          write (errmsg, '(a)') &
            'Simulation TIMING input file "'//trim(tdis6)// &
            '" does not exist.'
          call store_error(errmsg)
          call store_error_filename(simfile)
        end if
      case default
      end select
    end if
  end subroutine load_simtdis

  function remote_model_ndim(mtype, mfname) result(ncelldim)
    use SourceCommonModule, only: package_source_type
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: openfile, getunit
    use BlockParserModule, only: BlockParserType
    character(len=*), intent(in) :: mtype
    character(len=*), intent(in) :: mfname
    integer(I4B) :: ncelldim
    character(len=LENPACKAGENAME) :: source_type
    type(BlockParserType) :: parser
    integer(I4B) :: ierr, inunit
    logical(LGP) :: isfound, endOfBlock
    character(len=LINELENGTH) :: ptype

    ! initialize
    ncelldim = 0

    ! set source type
    source_type = package_source_type(mfname)
    select case (source_type)
    case ('MF6FILE')
      ! open name file
      inunit = getunit()
      call openfile(inunit, 0, trim(adjustl(mfname)), mtype, &
                    'FORMATTED', 'SEQUENTIAL', 'OLD')
      ! initialize parser
      call parser%Initialize(inunit, 0)
      ! get options block
      call parser%GetBlock('OPTIONS', isfound, ierr, &
                           supportOpenClose=.true., blockRequired=.false.)
      ! iterate through options
      if (isfound) then
        do
          call parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
        end do
      end if
      ! get packages block
      call parser%GetBlock('PACKAGES', isfound, ierr, &
                           supportOpenClose=.true., blockRequired=.true.)
      if (isfound) then
        ! read through packages
        do
          call parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          call parser%GetStringCaps(ptype)

          select case (ptype)
          case ('DIS6')
            ncelldim = 3
            exit
          case ('DIS2D6')
            ncelldim = 2
            exit
          case ('DISV6')
            ncelldim = 2
            exit
          case ('DISU6')
            ncelldim = 1
            exit
          case default
            write (errmsg, '(a)') &
              'Unknown discretization type "'//trim(ptype)// &
              '" not currently supported.'
            call store_error(errmsg, .true.)
          end select
        end do
      end if

      call parser%clear()
    case default
    end select
  end function remote_model_ndim

  !> @brief create model exports list
  !<
  subroutine export_cr()
    use ModelExportModule, only: modelexports_create, nc_export_active
#if defined(__WITH_NETCDF__)
    use NCExportCreateModule, only: nc_export_create
#endif
    call modelexports_create(iout)
    ! are netcdf exports elected
    if (nc_export_active()) then
#if defined(__WITH_NETCDF__)
      call nc_export_create()
#else
      write (errmsg, '(a)') &
        'Model namefile NETCDF_STUCTURED or NETCDF_MESH2D option configured &
        &but NetCDF libraries are not available.'
      call store_error(errmsg, .true.)
#endif
    end if
  end subroutine export_cr

  !> @brief model exports post prepare step actions
  !<
  subroutine export_post_prepare()
    use ModelExportModule, only: modelexports_post_prepare
    call modelexports_post_prepare()
  end subroutine export_post_prepare

  !> @brief model exports post step actions
  !<
  subroutine export_post_step()
    use ModelExportModule, only: modelexports_post_step
    call modelexports_post_step()
  end subroutine export_post_step

  !> @brief deallocate model export objects and list
  !<
  subroutine export_da()
    use ModelExportModule, only: modelexports_destroy
    call modelexports_destroy()
  end subroutine export_da

  !> @brief close an open netcdf file
  !<
  subroutine nc_close(ncid, nc_fname)
#if defined(__WITH_NETCDF__)
    use NetCDFCommonModule, only: nc_fclose
#endif
    integer(I4B), intent(in) :: ncid
    character(len=*), intent(in) :: nc_fname
    if (ncid > 0) then
#if defined(__WITH_NETCDF__)
      call nc_fclose(ncid, nc_fname)
#endif
    end if
  end subroutine nc_close

  !> @brief create model netcdf context
  !<
  function netcdf_context(modeltype, component_type, modelname, &
                          modelfname, iout) result(nc_vars)
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    use SourceCommonModule, only: filein_fname
#if defined(__WITH_NETCDF__)
    use NCContextBuildModule, only: open_ncfile, create_netcdf_context
#endif
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modelfname
    integer(I4B), intent(in) :: iout
    type(NCFileVarsType), pointer :: nc_vars
    character(len=LENMEMPATH) :: input_mempath
    character(len=LINELENGTH) :: nc_fname
    integer(I4B) :: ncid

    ! set input memory path
    input_mempath = create_mem_path(modelname, 'NAM', idm_context)

    ! allocate context object
    allocate (nc_vars)

    ! check if optional netcdf input file was provided
    if (filein_fname(nc_fname, 'NETCDF_FNAME', input_mempath, modelfname)) then
#if defined(__WITH_NETCDF__)
      ! open nc input file
      ncid = open_ncfile(nc_fname, iout)
      ! read the file and build the context
      call create_netcdf_context(modeltype, modelname, modelfname, &
                                 nc_vars, nc_fname, ncid, iout)
#else
      write (errmsg, '(a)') &
        'Cannot load model packages. NetCDF &
        &keyword specified in input file but &
        &NetCDF libraries are not available.'
      call store_error(errmsg)
      call store_error_filename(modelfname)
#endif
    else
      ncid = 0
      call nc_vars%init(modelname, '', ncid, '')
    end if
  end function netcdf_context

end module SourceLoadModule
