!> @brief This module contains the NCFileVarsModule
!!
!! These data structures organize package input information
!! associated with a single model netcdf input file.
!!
!<
module NCFileVarsModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENMODELNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use ListModule, only: ListType

  implicit none
  private
  public :: NCFileVarsType
  public :: NCPackageVarsType

  !> @brief Type describing input variables for a package in NetCDF file
  !<
  type :: NCPackageVarsType
    character(len=LENMODELNAME) :: modelname !< name of model
    type(ListType) :: nc_vars
    character(len=LINELENGTH), pointer :: grid => null() !< grid type
    character(len=LINELENGTH), pointer :: nc_fname => null() !< netcdf filename
    integer(I4B), pointer :: ncid => null() !< netcdf file handle
  contains
    procedure :: init => ncvars_init
    procedure :: varid => ncvars_varid
    procedure :: destroy => ncvars_destroy
  end type NCPackageVarsType

  !> @brief Type which describes a modflow input variable in a netcdf file
  !<
  type :: NCFileMf6VarType
    character(LINELENGTH) :: pkgname !< package name
    character(LINELENGTH) :: tagname !< tag name
    integer(I4B) :: layer !< variable layer
    integer(I4B) :: iaux !< variable aux index
    integer(I4B) :: varid !< NC file variable id
  contains
  end type NCFileMf6VarType

  !> @brief Type describing modflow6 input variables in model NetCDF file
  !<
  type :: NCFileVarsType
    type(ListType) :: mf6invar !< list of modflow 6 input variables in netcdf file
    character(len=LINELENGTH), pointer :: grid => null() !< grid type
    character(len=LINELENGTH), pointer :: nc_fname => null() !< netcdf filename
    integer(I4B), pointer :: ncid => null() !< netcdf file handle
  contains
    procedure :: init => fv_init
    procedure :: add => fv_add
    procedure :: destroy => fv_destroy
    procedure :: create_varlists
  end type NCFileVarsType

contains

  !> @brief create netcdf package variable lists
  !<
  subroutine ncvars_init(this, modelname)
    class(NCPackageVarsType) :: this
    character(len=*), intent(in) :: modelname
    ! set modelname
    this%modelname = modelname
  end subroutine ncvars_init

  !> @brief return a netcdf variable id for a package tagname
  !<
  function ncvars_varid(this, tagname, layer, iaux) result(varid)
    class(NCPackageVarsType) :: this
    character(len=*), intent(in) :: tagname
    integer(I4B), optional :: layer
    integer(I4B), optional :: iaux
    integer(I4B) :: varid
    integer(I4B) :: n, l, p, a
    class(NCFileMf6VarType), pointer :: nc_var

    ! initialize
    varid = -1
    l = -1
    p = -1
    a = -1

    ! set search layer if provided
    if (present(layer)) then
      l = layer
    end if

    ! set search iaux if provided
    if (present(iaux)) then
      a = iaux
    end if

    do n = 1, this%nc_vars%Count()
      nc_var => ncvar_get(this%nc_vars, n)
      if (nc_var%tagname == tagname .and. &
          nc_var%layer == l .and. &
          nc_var%iaux == a) then
        varid = nc_var%varid
      end if
    end do

    ! set error and exit if variable not in NetCDF input
    if (varid == -1) then
      if (this%nc_fname /= '') then
        write (errmsg, '(a)') &
          'NetCDF variable not found, tagname="'//trim(tagname)//'"'
        if (present(layer)) then
          write (errmsg, '(a,i0)') trim(errmsg)//', layer=', layer
        end if
        if (present(iaux)) then
          write (errmsg, '(a,i0)') trim(errmsg)//', iaux=', iaux
        end if
        write (errmsg, '(a)') trim(errmsg)//'.'
        call store_error(errmsg)
        call store_error_filename(this%nc_fname)
      else
        write (errmsg, '(a)') &
          'NetCDF variable not found, tagname="'//trim(tagname)// &
          '". NetCDF input not provided for model "'//trim(this%modelname)//'".'
        call store_error(errmsg, .true.)
      end if
    end if
  end function ncvars_varid

  !> @brief destroy netcdf package variable lists
  !<
  subroutine ncvars_destroy(this)
    class(NCPackageVarsType) :: this
    class(NCFileMf6VarType), pointer :: nc_var
    integer(I4B) :: n
    ! deallocate allocated memory
    do n = 1, this%nc_vars%Count()
      nc_var => ncvar_get(this%nc_vars, n)
      deallocate (nc_var)
      nullify (nc_var)
    end do
    call this%nc_vars%Clear()
  end subroutine ncvars_destroy

  !> @brief initialize netcdf model variable description type
  !<
  subroutine fv_init(this, modelname, nc_fname, ncid, grid)
    use ConstantsModule, only: LENMEMPATH
    use MemoryManagerModule, only: mem_allocate
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    class(NCFileVarsType) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: ncid
    character(len=*), intent(in) :: grid
    character(len=LENMEMPATH) :: mempath
    integer(I4B) :: ilen

    ! set mempath
    mempath = create_mem_path(component=modelname, &
                              context=idm_context)
    ! initialize strlen
    ilen = LINELENGTH

    ! allocate managed memory
    call mem_allocate(this%grid, ilen, 'NETCDF_GRID', mempath)
    call mem_allocate(this%nc_fname, ilen, 'NETCDF_FNAME', mempath)
    call mem_allocate(this%ncid, 'NCID', mempath)

    ! set
    this%grid = grid
    this%nc_fname = nc_fname
    this%ncid = ncid
  end subroutine fv_init

  !> @brief add netcdf modflow6 input variable to list
  !<
  subroutine fv_add(this, pkgname, tagname, layer, iaux, varid)
    use ArrayHandlersModule, only: expandarray
    class(NCFileVarsType) :: this
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    integer(I4B), intent(in) :: layer
    integer(I4B), intent(in) :: iaux
    integer(I4B), intent(in) :: varid
    class(NCFileMf6VarType), pointer :: invar
    class(*), pointer :: obj
    ! add mf6 variable to file list
    allocate (invar)
    invar%pkgname = pkgname
    invar%tagname = tagname
    invar%layer = layer
    invar%iaux = iaux
    invar%varid = varid
    obj => invar
    call this%mf6invar%Add(obj)
  end subroutine fv_add

  !> @brief destroy netcdf model variable description type
  !<
  subroutine fv_destroy(this)
    class(NCFileVarsType) :: this
    class(NCFileMf6VarType), pointer :: invar
    integer(I4B) :: n
    do n = 1, this%mf6invar%Count()
      invar => ncvar_get(this%mf6invar, n)
      deallocate (invar)
      nullify (invar)
    end do
    call this%mf6invar%Clear()
  end subroutine fv_destroy

  !> @brief create list of variables that correspond to a package
  !<
  subroutine create_varlists(this, modelname, pkgname, nc_vars)
    class(NCFileVarsType) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: pkgname
    type(NCPackageVarsType), pointer, intent(inout) :: nc_vars
    integer(I4B) :: n
    class(NCFileMf6VarType), pointer :: invar, nc_var
    class(*), pointer :: obj

    do n = 1, this%mf6invar%count()
      invar => ncvar_get(this%mf6invar, n)
      if (invar%pkgname == pkgname) then
        ! create package variable description
        allocate (nc_var)
        nc_var%pkgname = invar%pkgname
        nc_var%tagname = invar%tagname
        nc_var%layer = invar%layer
        nc_var%iaux = invar%iaux
        nc_var%varid = invar%varid
        obj => nc_var
        call nc_vars%nc_vars%Add(obj)
      end if
    end do

    ! set modelname
    nc_vars%modelname = modelname

    ! set file attribute pointers
    nc_vars%ncid => this%ncid
    nc_vars%nc_fname => this%nc_fname
    nc_vars%grid => this%grid
  end subroutine create_varlists

  !> @brief get modflow6 input variable description at position idx
  !<
  function ncvar_get(nc_vars, idx) result(res)
    type(ListType) :: nc_vars
    integer(I4B), intent(in) :: idx !< package number
    class(NCFileMf6VarType), pointer :: res
    class(*), pointer :: obj

    ! initialize res
    res => null()

    ! get the package from the list
    obj => nc_vars%GetItem(idx)
    if (associated(obj)) then
      select type (obj)
      class is (NCFileMf6VarType)
        res => obj
      end select
    end if
  end function ncvar_get

end module NCFileVarsModule
