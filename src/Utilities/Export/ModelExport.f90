!> @brief This module contains the ModelExportModule
!!
!! This modeule defines the local list of model exports.
!! It is not dependent on netcdf or other external
!! libraries.
!!
!<
module ModelExportModule

  use KindModule, only: DP, I4B, LGP
  use SimModule, only: store_error, store_error_filename
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH, LENMODELNAME, LENCOMPONENTNAME, &
                             LENMEMPATH
  use ListModule, only: ListType
  use NCModelExportModule, only: NCBaseModelExportType
  use InputLoadTypeModule, only: ModelDynamicPkgsType

  implicit none
  private
  public :: modelexports_create
  public :: modelexports_post_prepare
  public :: modelexports_post_step
  public :: modelexports_destroy
  public :: nc_export_active
  public :: export_models
  public :: get_export_model
  public :: ExportModelType

  type(ListType) :: export_models

  !> @brief export model type
  !!
  !!  This is a container variable which holds
  !!  model export objects.
  !!
  !<
  type :: ExportModelType
    type(ModelDynamicPkgsType), pointer :: loaders => null()
    character(len=LENMODELNAME) :: modelname !< name of model
    character(len=LENCOMPONENTNAME) :: modeltype !< type of model
    character(len=LINELENGTH) :: modelfname !< name of model input file
    character(len=LINELENGTH) :: nc_fname !< name of netcdf export file
    class(NCBaseModelExportType), pointer :: nc_export => null() !< netcdf export object pointer
    integer(I4B) :: nctype !< type of netcdf export
    integer(I4B) :: disenum !< type of discretization
    integer(I4B) :: iout !< lst file descriptor
  contains
    procedure :: init
    procedure :: post_prepare
    procedure :: post_step
    procedure :: destroy
  end type ExportModelType

contains

  !> @brief is netcdf export configured for any model
  !!
  function nc_export_active() result(active)
    use NCModelExportModule, only: NETCDF_UNDEF
    logical(LGP) :: active
    integer(I4B) :: n
    type(ExportModelType), pointer :: export_model
    active = .false.
    do n = 1, export_models%Count()
      export_model => get_export_model(n)
      if (export_model%nctype /= NETCDF_UNDEF) then
        active = .true.
        exit
      end if
    end do
  end function nc_export_active

  !> @brief create export container variable for all local models
  !!
  subroutine modelexports_create(iout)
    use InputLoadTypeModule, only: model_inputs
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryHelperModule, only: create_mem_path
    use InputLoadTypeModule, only: GetDynamicModelFromList
    use SimVariablesModule, only: idm_context
    use NCModelExportModule, only: NETCDF_MESH2D, NETCDF_STRUCTURED
    use SourceCommonModule, only: file_ext
    integer(I4B), intent(in) :: iout
    type(ModelDynamicPkgsType), pointer :: model_input
    type(ExportModelType), pointer :: export_model
    character(len=LENMEMPATH) :: modelnam_mempath, model_mempath, ext
    integer(I4B), pointer :: disenum
    integer(I4B) :: n
    logical(LGP) :: found

    do n = 1, model_inputs%Count()
      ! allocate and initialize
      allocate (export_model)

      ! set pointer to dynamic input model instance
      model_input => GetDynamicModelFromList(model_inputs, n)

      ! set input mempaths
      modelnam_mempath = &
        create_mem_path(component=model_input%modelname, &
                        subcomponent='NAM', context=idm_context)
      model_mempath = create_mem_path(component=model_input%modelname, &
                                      context=idm_context)
      ! set pointer to dis enum type
      call mem_setptr(disenum, 'DISENUM', model_mempath)

      ! initialize model
      call export_model%init(model_input, disenum, iout)

      ! update NetCDF fileout name if provided
      call mem_set_value(export_model%nc_fname, 'NCMESH2DFILE', &
                         modelnam_mempath, found)
      if (found) then
        export_model%nctype = NETCDF_MESH2D
      else
        call mem_set_value(export_model%nc_fname, 'NCSTRUCTFILE', &
                           modelnam_mempath, found)
        if (found) then
          export_model%nctype = NETCDF_STRUCTURED
        end if
      end if

      if (found) then
        ext = file_ext(export_model%nc_fname)
        if (ext /= 'nc') then
          errmsg = 'NetCDF output file name must use ".nc" extension. '// &
                   'Filename="'//trim(export_model%nc_fname)//'".'
          call store_error(errmsg)
          call store_error_filename(export_model%modelfname)
        end if
      end if

      ! add model to list
      call add_export_model(export_model)
    end do
  end subroutine modelexports_create

  !> @brief export model list post prepare step
  !!
  subroutine modelexports_post_prepare()
    class(*), pointer :: obj
    class(ExportModelType), pointer :: export_model
    integer(I4B) :: n
    do n = 1, export_models%Count()
      obj => export_models%GetItem(n)
      if (associated(obj)) then
        select type (obj)
        class is (ExportModelType)
          export_model => obj
          call export_model%post_prepare()
        end select
      end if
    end do
  end subroutine modelexports_post_prepare

  !> @brief export model list post step
  !!
  subroutine modelexports_post_step()
    class(*), pointer :: obj
    class(ExportModelType), pointer :: export_model
    integer(I4B) :: n
    do n = 1, export_models%Count()
      obj => export_models%GetItem(n)
      if (associated(obj)) then
        select type (obj)
        class is (ExportModelType)
          export_model => obj
          call export_model%post_step()
        end select
      end if
    end do
  end subroutine modelexports_post_step

  !> @brief destroy export model list
  !!
  subroutine modelexports_destroy()
    class(*), pointer :: obj
    class(ExportModelType), pointer :: export_model
    integer(I4B) :: n
    do n = 1, export_models%Count()
      obj => export_models%GetItem(n)
      if (associated(obj)) then
        select type (obj)
        class is (ExportModelType)
          export_model => obj
          call export_model%destroy()
          deallocate (export_model)
          nullify (export_model)
        end select
      end if
    end do
    !
    call export_models%clear()
  end subroutine modelexports_destroy

  !> @brief initialize model export container variable
  !!
  !<
  subroutine init(this, loaders, disenum, iout)
    use NCModelExportModule, only: NETCDF_UNDEF
    class(ExportModelType), intent(inout) :: this
    type(ModelDynamicPkgsType), pointer, intent(in) :: loaders
    integer(I4B), intent(in) :: disenum
    integer(I4B), intent(in) :: iout
    this%loaders => loaders
    this%modelname = loaders%modelname
    this%modeltype = loaders%modeltype
    this%modelfname = loaders%modelfname
    this%nc_fname = ''
    this%nctype = NETCDF_UNDEF
    this%disenum = disenum
    this%iout = iout
    nullify (this%nc_export)
  end subroutine init

  !> @brief model export container post prepare step actions
  !!
  !<
  subroutine post_prepare(this)
    class(ExportModelType), intent(inout) :: this
    if (associated(this%nc_export)) then
      call this%nc_export%export_input()
    end if
  end subroutine post_prepare

  !> @brief model export container post step actions
  !!
  !<
  subroutine post_step(this)
    class(ExportModelType), intent(inout) :: this
    if (associated(this%nc_export)) then
      call this%nc_export%step()
    end if
  end subroutine post_step

  !> @brief destroy model export container
  !!
  !<
  subroutine destroy(this)
    class(ExportModelType), intent(inout) :: this
    if (associated(this%nc_export)) then
      call this%nc_export%destroy()
      deallocate (this%nc_export)
      nullify (this%nc_export)
    end if
  end subroutine destroy

  !> @brief add model export object to list
  !!
  !<
  subroutine add_export_model(export_model)
    type(ExportModelType), pointer, intent(inout) :: export_model
    class(*), pointer :: obj
    obj => export_model
    call export_models%Add(obj)
  end subroutine add_export_model

  !> @brief get model export object by index
  !!
  !<
  function get_export_model(idx) result(res)
    integer(I4B), intent(in) :: idx !< package number
    class(ExportModelType), pointer :: res
    class(*), pointer :: obj
    ! initialize res
    nullify (res)
    ! get the object from the list
    obj => export_models%GetItem(idx)
    if (associated(obj)) then
      select type (obj)
      class is (ExportModelType)
        res => obj
      end select
    end if
  end function get_export_model

end module ModelExportModule
