!> @brief This module contains the ModelPackageInputsModule
!!
!! This module contains the high-level routines for assembling
!! model package information and loading to the input context
!!
!<
module ModelPackageInputsModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENMODELNAME, LENFTYPE, &
                             LENPACKAGETYPE, LENPACKAGENAME, LENCOMPONENTNAME
  use SimModule, only: store_error, count_errors, store_error_filename
  use SimVariablesModule, only: iout
  use ArrayHandlersModule, only: expandarray
  use CharacterStringModule, only: CharacterStringType

  implicit none
  private
  public :: ModelPackageInputsType
  public :: supported_model

  !> @brief derived type for loadable package type
  !!
  !!  This derived type is used to store package instance
  !!  descriptions for a supported package type.
  !!
  !<
  type :: LoadablePackageType
    ! package type, e.g. 'DIS6' or 'CHD6'
    character(len=LENPACKAGETYPE) :: pkgtype
    ! component type, e.g. 'DIS' or 'CHD'
    character(len=LENCOMPONENTNAME) :: subcomponent_type
    ! package instance attribute arrays
    character(len=LINELENGTH), dimension(:), allocatable :: filenames
    character(len=LENPACKAGENAME), dimension(:), allocatable :: pkgnames
    character(len=LENMEMPATH), dimension(:), allocatable :: mempaths
    integer(I4B), dimension(:), allocatable :: inunits
    ! number of package instances
    integer(I4B) :: pnum
  contains
    procedure :: create => pkgtype_create
    procedure :: add => pkgtype_add
    procedure :: destroy => pkgtype_destroy
  end type LoadablePackageType

  !> @brief derived type for model package inputs type
  !!
  !!  This derived type is used to define input package
  !!  descriptors for a model and load to managed memory.
  !!
  !<
  type :: ModelPackageInputsType
    ! model attributes
    character(len=LENPACKAGETYPE) :: modeltype ! model type, e.g. 'GWF6'
    character(len=LINELENGTH) :: modelfname
    character(len=LENMODELNAME) :: modelname
    ! component type
    character(len=LENCOMPONENTNAME) :: component_type ! e.g. 'GWF'
    ! mempaths
    character(len=LENMEMPATH) :: input_mempath
    character(len=LENMEMPATH) :: model_mempath
    ! pointers to created managed memory
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    integer(I4B), dimension(:), contiguous, &
      pointer :: inunits => null()
    ! loadable package type array
    type(LoadablePackageType), dimension(:), allocatable :: pkglist
    ! pkgtype definitions
    integer(I4B) :: niunit
    character(len=LENPACKAGETYPE), dimension(:), allocatable :: cunit
    ! out handle
    integer(I4B) :: iout
  contains
    procedure :: init => modelpkgs_init
    procedure :: memload => modelpkgs_memload
    procedure :: destroy => modelpkgs_destroy
    procedure, private :: create => modelpkgs_create
    procedure, private :: addpkgs => modelpkgs_addpkgs
    procedure, private :: add => modelpkgs_add
    procedure, private :: pkgcount => modelpkgs_pkgcount
  end type ModelPackageInputsType

contains

  !> @brief does model support multiple instances of this package type
  !<
  function multi_pkg_type(mtype_component, ptype_component, pkgtype) &
    result(multi_pkg)
    use IdmDfnSelectorModule, only: idm_integrated, idm_multi_package
    use ModelPackageInputModule, only: multi_package_type
    character(len=LENCOMPONENTNAME), intent(in) :: mtype_component
    character(len=LENCOMPONENTNAME), intent(in) :: ptype_component
    character(len=LENFTYPE), intent(in) :: pkgtype
    logical(LGP) :: multi_pkg
    multi_pkg = .false.
    if (idm_integrated(mtype_component, ptype_component)) then
      multi_pkg = idm_multi_package(mtype_component, ptype_component)
    else
      multi_pkg = multi_package_type(mtype_component, ptype_component, pkgtype)
    end if
  end function multi_pkg_type

  !> @brief is this a supported MODFLOW 6 model type
  !<
  function supported_model(ctype)
    use ModelPackageInputModule, only: NMODEL, MODFLOW6MODELS
    character(len=*), intent(in) :: ctype
    logical(LGP) :: supported_model
    integer(I4B) :: n
    supported_model = .false.
    do n = 1, NMODEL
      if (ctype == MODFLOW6MODELS(n)) then
        supported_model = .true.
        exit
      end if
    end do
  end function supported_model

  !> @brief create a new package type
  !<
  subroutine pkgtype_create(this, modeltype, modelname, pkgtype)
    use SourceCommonModule, only: idm_subcomponent_type
    class(LoadablePackageType) :: this
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: pkgtype

    ! initialize
    this%pkgtype = pkgtype
    this%subcomponent_type = idm_subcomponent_type(modeltype, pkgtype)
    this%pnum = 0

    ! allocate arrays
    allocate (this%filenames(0))
    allocate (this%pkgnames(0))
    allocate (this%mempaths(0))
    allocate (this%inunits(0))
  end subroutine pkgtype_create

  !> @brief add a new package instance to this package type
  !<
  subroutine pkgtype_add(this, modelname, mtype_component, filetype, &
                         filename, pkgname, iout)
    use MemoryManagerModule, only: mem_allocate
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    use IdmDfnSelectorModule, only: idm_integrated, idm_multi_package
    use SourceCommonModule, only: idm_subcomponent_name
    class(LoadablePackageType) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: mtype_component
    character(len=*), intent(in) :: filetype
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: pkgname
    integer(I4B), intent(in) :: iout
    character(len=LENPACKAGENAME) :: sc_name, pname
    character(len=LENMEMPATH) :: mempath
    character(len=LINELENGTH), pointer :: cstr

    ! reallocate
    call expandarray(this%filenames)
    call expandarray(this%pkgnames)
    call expandarray(this%inunits)
    call expandarray(this%mempaths)

    ! add new package instance
    this%pnum = this%pnum + 1
    this%filenames(this%pnum) = filename
    this%pkgnames(this%pnum) = pkgname
    this%inunits(this%pnum) = 0

    ! set pkgname if empty
    if (this%pkgnames(this%pnum) == '') then
      if (multi_pkg_type(mtype_component, &
                         this%subcomponent_type, &
                         filetype)) then
        write (pname, '(a,i0)') trim(this%subcomponent_type)//'-', this%pnum
      else
        write (pname, '(a,i0)') trim(this%subcomponent_type)
      end if
      this%pkgnames(this%pnum) = pname
    end if

    ! set up input context for model
    if (idm_integrated(mtype_component, this%subcomponent_type)) then
      ! set subcomponent name
      sc_name = idm_subcomponent_name(mtype_component, this%subcomponent_type, &
                                      this%pkgnames(this%pnum))
      ! create and store the mempath
      this%mempaths(this%pnum) = &
        create_mem_path(modelname, sc_name, idm_context)
      ! allocate and initialize filename for package
      mempath = create_mem_path(modelname, sc_name, idm_context)
      call mem_allocate(cstr, LINELENGTH, 'INPUT_FNAME', mempath)
      cstr = filename
    else
      ! set mempath empty
      this%mempaths(this%pnum) = ''
    end if
  end subroutine pkgtype_add

  !> @brief deallocate object
  !<
  subroutine pkgtype_destroy(this)
    class(LoadablePackageType) :: this
    ! deallocate dynamic arrays
    deallocate (this%filenames)
    deallocate (this%pkgnames)
    deallocate (this%inunits)
    deallocate (this%mempaths)
  end subroutine pkgtype_destroy

  !> @brief initialize model package inputs object
  !<
  subroutine modelpkgs_init(this, modeltype, modelfname, modelname, iout)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_allocate
    use SimVariablesModule, only: idm_context, simfile
    use SourceCommonModule, only: idm_component_type
    use ModelPackageInputModule, only: supported_model_packages
    class(ModelPackageInputsType) :: this
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: modelname
    integer(I4B), intent(in) :: iout

    ! initialize object
    this%modeltype = modeltype
    this%modelfname = modelfname
    this%modelname = modelname
    this%component_type = idm_component_type(modeltype)
    this%iout = iout

    ! verify user specified model type
    if (.not. supported_model(modeltype)) then
      ! -- error and exit for unsupported model type
      write (errmsg, '(3a)') 'Models block model type "', trim(modeltype), &
        '" is not valid.'
      call store_error(errmsg)
      call store_error_filename(simfile)
    end if

    ! allocate and set model supported package types
    call supported_model_packages(modeltype, this%cunit, this%niunit)

    ! set memory paths
    this%input_mempath = create_mem_path(this%modelname, 'NAM', idm_context)
    this%model_mempath = create_mem_path(component=this%modelname, &
                                         context=idm_context)
    ! allocate managed memory
    call mem_allocate(this%pkgtypes, LENPACKAGETYPE, 0, 'PKGTYPES', &
                      this%model_mempath)
    call mem_allocate(this%pkgnames, LENPACKAGENAME, 0, 'PKGNAMES', &
                      this%model_mempath)
    call mem_allocate(this%mempaths, LENMEMPATH, 0, 'MEMPATHS', &
                      this%model_mempath)
    call mem_allocate(this%inunits, 0, 'INUNITS', this%model_mempath)

    ! build descriptions of packages
    call this%addpkgs()
  end subroutine modelpkgs_init

  !> @brief create the package type list
  !<
  subroutine modelpkgs_create(this, ftypes)
    use SortModule, only: qsort
    class(ModelPackageInputsType) :: this
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: ftypes
    integer(I4B), dimension(:), allocatable :: cunit_idxs, indx
    character(len=LENPACKAGETYPE) :: ftype
    integer(I4B) :: n, m
    logical(LGP) :: found

    ! allocate
    allocate (cunit_idxs(0))

    ! identify input packages and check that each is supported
    do n = 1, size(ftypes)
      ! type from model nam file packages block
      ftype = ftypes(n)
      found = .false.

      ! search supported types for this filetype
      do m = 1, this%niunit
        if (this%cunit(m) == ftype) then
          ! set found
          found = .true.

          ! add to cunit list if first instance of this type
          if (any(cunit_idxs == m)) then
            ! no-op
          else
            call expandarray(cunit_idxs)
            cunit_idxs(size(cunit_idxs)) = m
          end if

          ! exit search
          exit
        end if
      end do

      ! set error if namfile pkg filetype is not supported
      if (.not. found) then
        write (errmsg, '(a,a,a,a,a)') 'Model package type not supported &
          &[model=', trim(this%modelname), ', type=', &
          trim(ftype), '].'
        call store_error(errmsg)
        call store_error_filename(this%modelfname)
      end if
    end do

    ! allocate the pkglist
    allocate (this%pkglist(size(cunit_idxs)))

    ! sort cunit indexes
    allocate (indx(size(cunit_idxs)))
    call qsort(indx, cunit_idxs)

    ! create sorted LoadablePackageType object list
    do n = 1, size(cunit_idxs)
      call this%pkglist(n)%create(this%modeltype, this%modelname, &
                                  this%cunit(cunit_idxs(n)))
    end do

    ! cleanup
    deallocate (cunit_idxs)
    deallocate (indx)
  end subroutine modelpkgs_create

  !> @brief add a model package instance to package type list
  !<
  subroutine modelpkgs_add(this, pkgtype, filename, pkgname)
    class(ModelPackageInputsType) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: pkgname
    type(LoadablePackageType) :: pkg
    integer(I4B) :: n
    ! locate index of pkgtype in pkglist
    do n = 1, size(this%pkglist)
      pkg = this%pkglist(n)
      if (pkg%pkgtype == pkgtype) then
        call this%pkglist(n)%add(this%modelname, this%component_type, &
                                 pkgtype, filename, pkgname, this%iout)
        exit
      end if
    end do
  end subroutine modelpkgs_add

  !> @brief build the type list with all model package instances
  !<
  subroutine modelpkgs_addpkgs(this)
    use MemoryManagerModule, only: mem_setptr
    use SourceCommonModule, only: inlen_check
    class(ModelPackageInputsType) :: this
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: ftypes !< file types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: fnames !< file names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pnames !< package names
    character(len=LINELENGTH) :: ftype, fname, pname
    integer(I4B) :: n

    ! set pointers to input context model package attribute arrays
    call mem_setptr(ftypes, 'FTYPE', this%input_mempath)
    call mem_setptr(fnames, 'FNAME', this%input_mempath)
    call mem_setptr(pnames, 'PNAME', this%input_mempath)

    ! create the package list
    call this%create(ftypes)

    ! load model packages
    do n = 1, size(ftypes)
      ! attributes for this package
      ftype = ftypes(n)
      fname = fnames(n)
      call inlen_check(pnames(n), pname, LENPACKAGENAME, 'PACKAGENAME')

      ! add this instance to package list
      call this%add(ftype, fname, pname)
    end do

    ! terminate if errors were detected
    if (count_errors() > 0) then
      call store_error_filename(this%modelfname)
    end if
  end subroutine modelpkgs_addpkgs

  !> @brief get package instance count and verify base or multi of each
  !<
  function modelpkgs_pkgcount(this) result(pnum)
    class(ModelPackageInputsType) :: this
    integer(I4B) :: pnum
    integer(I4B) :: n

    ! initialize
    pnum = 0

    ! count model package instances
    do n = 1, size(this%pkglist)
      if (multi_pkg_type(this%component_type, &
                         this%pkglist(n)%subcomponent_type, &
                         this%pkglist(n)%pkgtype)) then
        ! multiple instances ok
      else
        ! set error for unexpected extra packages
        if (this%pkglist(n)%pnum > 1) then
          write (errmsg, '(a,a,a,a,a)') &
            'Multiple instances specified for model base package type &
            &[model=', trim(this%modelname), ', type=', &
            trim(this%pkglist(n)%pkgtype), '].'
          call store_error(errmsg)
          call store_error_filename(this%modelfname)
        end if
      end if

      ! add to package count
      pnum = pnum + this%pkglist(n)%pnum
    end do
  end function modelpkgs_pkgcount

  !> @brief load package descriptors to managed memory
  !<
  subroutine modelpkgs_memload(this)
    use MemoryManagerModule, only: mem_reallocate
    class(ModelPackageInputsType) :: this
    integer(I4B) :: n, m, idx
    integer(I4B) :: pnum

    ! initialize load index
    idx = 0

    ! set total number of package instances
    pnum = this%pkgcount()

    ! reallocate model input package attribute arrays
    call mem_reallocate(this%pkgtypes, LENPACKAGETYPE, pnum, 'PKGTYPES', &
                        this%model_mempath)
    call mem_reallocate(this%pkgnames, LENPACKAGENAME, pnum, 'PKGNAMES', &
                        this%model_mempath)
    call mem_reallocate(this%mempaths, LENMEMPATH, pnum, 'MEMPATHS', &
                        this%model_mempath)
    call mem_reallocate(this%inunits, pnum, 'INUNITS', this%model_mempath)

    ! load pkinfo
    do n = 1, size(this%pkglist)
      do m = 1, this%pkglist(n)%pnum
        ! increment index
        idx = idx + 1
        ! package type like 'CHD6'
        this%pkgtypes(idx) = trim(this%pkglist(n)%pkgtype)
        ! package name like 'CHD-2'
        this%pkgnames(idx) = trim(this%pkglist(n)%pkgnames(m))
        ! memory path like '__INPUT__/MYMODEL/CHD-2'
        this%mempaths(idx) = trim(this%pkglist(n)%mempaths(m))
        ! input file unit number
        this%inunits(idx) = this%pkglist(n)%inunits(m)
      end do
    end do
  end subroutine modelpkgs_memload

  !> @brief deallocate object
  !<
  subroutine modelpkgs_destroy(this)
    class(ModelPackageInputsType) :: this
    integer(I4B) :: n
    do n = 1, size(this%pkglist)
      call this%pkglist(n)%destroy()
    end do
    deallocate (this%pkglist)
    deallocate (this%cunit)
  end subroutine modelpkgs_destroy

end module ModelPackageInputsModule
