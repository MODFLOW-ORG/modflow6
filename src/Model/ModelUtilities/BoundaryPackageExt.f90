!> @brief This module contains the extended boundary package
!!
!! This module contains the extended boundary type that itself
!! should be extended by model boundary packages that have been
!! updated to source static and dynamic input data from the
!! input context.
!!
!<
module BndExtModule

  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: LENMEMPATH, LENBOUNDNAME, LENAUXNAME, LINELENGTH
  use ObsModule, only: obs_cr
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors, store_error_filename
  use BndModule, only: BndType
  use GeomUtilModule, only: get_node, get_ijk

  implicit none

  private
  public :: BndExtType

  !> @ brief BndExtType
  !!
  !!  Generic extended boundary package type.  This derived type can be
  !!  overridden to define concrete boundary package types that source
  !!  all input from the input context.
  !<
  type, extends(BndType) :: BndExtType
    ! -- characters
    ! -- scalars
    integer(I4B), pointer :: iper
    logical(LGP), pointer :: readarraygrid
    logical(LGP), pointer :: readarraylayer
    ! -- arrays
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid => null() !< input user cellid list
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist => null() !< input user nodelist
  contains
    procedure :: bnd_df => bndext_df
    procedure :: bnd_rp => bndext_rp
    procedure :: bnd_da => bndext_da
    procedure :: allocate_scalars => bndext_allocate_scalars
    procedure :: allocate_arrays => bndext_allocate_arrays
    procedure :: source_options
    procedure :: source_dimensions
    procedure :: log_options
    procedure :: nodelist_update
    procedure :: check_cellid
    procedure :: write_list
    procedure :: bound_value
  end type BndExtType

  !> @ brief BndExtFoundType
  !!
  !!  This type is used to simplify the tracking of common parameters
  !!  that are sourced from the input context.
  !<
  type BndExtFoundType
    logical :: naux = .false.
    logical :: ipakcb = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: boundnames = .false.
    logical :: auxmultname = .false.
    logical :: inewton = .false.
    logical :: auxiliary = .false.
    logical :: maxbound = .false.
  end type BndExtFoundType

contains

  !> @ brief Define boundary package options and dimensions
    !!
    !!  Define base boundary package options and dimensions for
    !!  a model boundary package.
    !!
  !<
  subroutine bndext_df(this, neq, dis)
    ! -- modules
    use BaseDisModule, only: DisBaseType
    use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType, &
                                            tasmanager_cr
    use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
    ! -- dummy variables
    class(BndExtType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(inout) :: neq !< number of equations
    class(DisBaseType), pointer :: dis !< discretization object
    !
    ! -- set pointer to dis object for the model
    this%dis => dis
    !
    ! -- Create time series managers
    ! -- Not in use by this type but BndType uses and deallocates
    call tsmanager_cr(this%TsManager, this%iout)
    call tasmanager_cr(this%TasManager, dis, this%name_model, this%iout)
    !
    ! -- create obs package
    call obs_cr(this%obs, this%inobspkg)
    !
    ! -- Write information to model list file
    write (this%iout, 1) this%filtyp, trim(adjustl(this%text)), this%input_mempath
1   format(1X, /1X, a, ' -- ', a, ' PACKAGE, VERSION 8, 2/22/2014', &
           ' INPUT READ FROM MEMPATH: ', a)
    !
    ! -- source options
    call this%source_options()
    !
    ! -- Define time series managers
    call this%tsmanager%tsmanager_df()
    call this%tasmanager%tasmanager_df()
    !
    ! -- source dimensions
    call this%source_dimensions()
    !
    ! -- update package moffset for packages that add rows
    if (this%npakeq > 0) then
      this%ioffset = neq - this%dis%nodes
    end if
    !
    ! -- update neq
    neq = neq + this%npakeq
    !
    ! -- Store information needed for observations
    if (this%bnd_obs_supported()) then
      call this%obs%obs_df(this%iout, this%packName, this%filtyp, this%dis)
      call this%bnd_df_obs()
    end if
  end subroutine bndext_df

  subroutine bndext_rp(this)
    ! -- modules
    use TdisModule, only: kper
    use MemoryManagerModule, only: mem_deallocate, mem_reallocate
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy variables
    class(BndExtType), intent(inout) :: this !< BndExtType object
    ! -- local variables
    logical(LGP) :: found
    integer(I4B) :: n, noder, nodeuser
    character(len=LINELENGTH) :: nodestr
    !
    if (this%iper /= kper) return
    !
    ! -- copy nbound from input context
    call mem_set_value(this%nbound, 'NBOUND', this%input_mempath, &
                       found)

    if (this%readarraygrid) then
      ! -- Set the nodelist
      do n = 1, this%nbound
        nodeuser = this%nodeulist(n)
        noder = this%dis%get_nodenumber(nodeuser, 1)
        if (noder >= 0) then
          this%nodelist(n) = noder
        else
          call this%dis%nodeu_to_string(n, nodestr)
          write (errmsg, *) &
            ' Cell is outside active grid domain: '// &
            trim(adjustl(nodestr))
          call store_error(errmsg)
        end if
      end do
      !
      ! -- exit if errors were found
      if (count_errors() > 0) then
        write (errmsg, *) count_errors(), ' errors encountered.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
    else
      !
      ! -- convert cellids to node numbers
      call this%nodelist_update()
      !
      ! -- update boundname string list
      if (this%inamedbound /= 0) then
        do n = 1, size(this%boundname_cst)
          this%boundname(n) = this%boundname_cst(n)
        end do
      end if
    end if
  end subroutine bndext_rp

  !> @ brief Deallocate package memory
  !<
  subroutine bndext_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate, mem_setptr
    ! -- dummy variables
    class(BndExtType) :: this !< BndExtType object
    !
    ! -- deallocate checkin paths
    call mem_deallocate(this%cellid, 'CELLID', this%memoryPath)
    call mem_deallocate(this%nodeulist, 'NODEULIST', this%memoryPath)
    call mem_deallocate(this%boundname_cst, 'BOUNDNAME_IDM', this%memoryPath)
    call mem_deallocate(this%auxvar, 'AUXVAR_IDM', this%memoryPath)
    !
    ! -- reassign pointers for base class _da
    call mem_setptr(this%boundname_cst, 'BOUNDNAME_CST', this%memoryPath)
    call mem_setptr(this%auxvar, 'AUXVAR', this%memoryPath)
    !
    ! -- scalars
    deallocate (this%readarraygrid)
    deallocate (this%readarraylayer)
    nullify (this%readarraygrid)
    nullify (this%readarraylayer)
    nullify (this%iper)
    !
    ! -- deallocate
    call this%BndType%bnd_da()
  end subroutine bndext_da

  !> @ brief Allocate package scalars
    !!
    !!  Allocate and initialize base boundary package scalars. This method
    !!  only needs to be overridden if additional scalars are defined
    !!  for a specific package.
    !!
  !<
  subroutine bndext_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    ! -- dummy variables
    class(BndExtType) :: this !< BndExtType object
    ! -- local variables
    character(len=LENMEMPATH) :: input_mempath
    logical(LGP) :: found
    !
    ! -- set memory path
    input_mempath = create_mem_path(this%name_model, this%packName, idm_context)
    !
    ! -- allocate base BndType scalars
    call this%BndType%allocate_scalars()
    !
    ! -- set IPER pointer
    call mem_setptr(this%iper, 'IPER', input_mempath)

    ! -- allocate internal scalars
    allocate (this%readarraygrid)
    allocate (this%readarraylayer)

    ! -- initialize internal scalars
    this%readarraygrid = .false.
    this%readarraylayer = .false.

    ! -- update internal scalars based on user input
    call mem_set_value(this%readarraygrid, 'READARRAYGRID', input_mempath, found)
    call mem_set_value(this%readarraylayer, 'READARRAYLAYER', &
                       input_mempath, found)

    ! -- no packages currently use READARRAYLAYER
    if (this%readarraylayer) then
      write (errmsg, '(a)') 'READARRAYLAYER is not currently supported.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
  end subroutine bndext_allocate_scalars

  !> @ brief Allocate package arrays
    !!
    !!  Allocate and initialize base boundary package arrays. This method
    !!  only needs to be overridden if additional arrays are defined
    !!  for a specific package.
    !!
  !<
  subroutine bndext_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate, mem_setptr, mem_checkin
    ! -- dummy variables
    class(BndExtType) :: this !< BndExtType object
    ! -- local variables
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist !< package nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar !< package aux variable array
    !
    ! -- allocate base BndType arrays
    call this%BndType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set input context pointers
    call mem_setptr(this%cellid, 'CELLID', this%input_mempath)
    call mem_setptr(this%nodeulist, 'NODEULIST', this%input_mempath)
    call mem_setptr(this%boundname_cst, 'BOUNDNAME', this%input_mempath)
    !
    ! -- checkin input context pointers
    call mem_checkin(this%cellid, 'CELLID', this%memoryPath, &
                     'CELLID', this%input_mempath)
    call mem_checkin(this%nodeulist, 'NODEULIST', this%memoryPath, &
                     'NODEULIST', this%input_mempath)
    call mem_checkin(this%boundname_cst, LENBOUNDNAME, 'BOUNDNAME_IDM', &
                     this%memoryPath, 'BOUNDNAME', this%input_mempath)
    !
    if (present(auxvar)) then
      ! no-op
    else
      ! -- set auxvar input context pointer
      call mem_setptr(this%auxvar, 'AUXVAR', this%input_mempath)
      !
      ! -- checkin auxvar input context pointer
      call mem_checkin(this%auxvar, 'AUXVAR_IDM', this%memoryPath, &
                       'AUXVAR', this%input_mempath)
    end if
  end subroutine bndext_allocate_arrays

  !> @ brief Source package options from input context
  !<
  subroutine source_options(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate, mem_setptr !, get_isize
    use MemoryManagerExtModule, only: mem_set_value
    use InputOutputModule, only: GetUnit, openfile
    use CharacterStringModule, only: CharacterStringType
    use SourceCommonModule, only: filein_fname
    ! -- dummy variables
    class(BndExtType), intent(inout) :: this !< BndExtType object
    ! -- local variables
    type(BndExtFoundType) :: found
    character(len=LENAUXNAME) :: sfacauxname
    integer(I4B) :: n
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%naux, 'NAUX', this%input_mempath, found%naux)
    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)
    call mem_set_value(this%iprpak, 'IPRPAK', this%input_mempath, found%iprpak)
    call mem_set_value(this%iprflow, 'IPRFLOW', this%input_mempath, found%iprflow)
    call mem_set_value(this%inamedbound, 'BOUNDNAMES', this%input_mempath, &
                       found%boundnames)
    call mem_set_value(sfacauxname, 'AUXMULTNAME', this%input_mempath, &
                       found%auxmultname)
    call mem_set_value(this%inewton, 'INEWTON', this%input_mempath, found%inewton)
    !
    ! -- log found options
    call this%log_options(found, sfacauxname)
    !
    ! -- reallocate aux arrays if aux variables provided
    if (found%naux .and. this%naux > 0) then
      call mem_reallocate(this%auxname, LENAUXNAME, this%naux, &
                          'AUXNAME', this%memoryPath)
      call mem_reallocate(this%auxname_cst, LENAUXNAME, this%naux, &
                          'AUXNAME_CST', this%memoryPath)
      call mem_set_value(this%auxname_cst, 'AUXILIARY', this%input_mempath, &
                         found%auxiliary)
      !
      do n = 1, this%naux
        this%auxname(n) = this%auxname_cst(n)
      end do
    end if
    !
    ! -- save flows option active
    if (found%ipakcb) this%ipakcb = -1
    !
    ! -- auxmultname provided
    if (found%auxmultname) this%iauxmultcol = -1
    !
    !
    ! -- enforce 0 or 1 OBS6_FILENAME entries in option block
    if (filein_fname(this%obs%inputFilename, 'OBS6_FILENAME', &
                     this%input_mempath, this%input_fname)) then
      this%obs%active = .true.
      this%obs%inUnitObs = GetUnit()
      call openfile(this%obs%inUnitObs, this%iout, this%obs%inputFilename, 'OBS')
    end if
    !
    ! -- no newton specified
    if (found%inewton) this%inewton = 0
    !
    ! -- AUXMULTNAME was specified, so find column of auxvar that will be multiplier
    if (this%iauxmultcol < 0) then
      !
      ! -- Error if no aux variable specified
      if (this%naux == 0) then
        write (errmsg, '(a,2(1x,a))') &
          'AUXMULTNAME was specified as', trim(adjustl(sfacauxname)), &
          'but no AUX variables specified.'
        call store_error(errmsg)
      end if
      !
      ! -- Assign mult column
      this%iauxmultcol = 0
      do n = 1, this%naux
        if (sfacauxname == this%auxname(n)) then
          this%iauxmultcol = n
          exit
        end if
      end do
      !
      ! -- Error if aux variable cannot be found
      if (this%iauxmultcol == 0) then
        write (errmsg, '(a,2(1x,a))') &
          'AUXMULTNAME was specified as', trim(adjustl(sfacauxname)), &
          'but no AUX variable found with this name.'
        call store_error(errmsg)
      end if
    end if
    !
    ! -- terminate if errors were detected
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine source_options

  !> @ brief Log package options
  !<
  subroutine log_options(this, found, sfacauxname)
    ! -- modules
    ! -- dummy variables
    class(BndExtType), intent(inout) :: this !< BndExtType object
    type(BndExtFoundType), intent(in) :: found
    character(len=*), intent(in) :: sfacauxname
    ! -- local variables
    ! -- format
    character(len=*), parameter :: fmtflow = &
      &"(4x, 'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL')"
    character(len=*), parameter :: fmttas = &
      &"(4x, 'TIME-ARRAY SERIES DATA WILL BE READ FROM FILE: ', a)"
    character(len=*), parameter :: fmtts = &
      &"(4x, 'TIME-SERIES DATA WILL BE READ FROM FILE: ', a)"
    character(len=*), parameter :: fmtnme = &
      &"(a, i0, a)"
    !
    ! -- log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' BASE OPTIONS'
    !
    if (found%ipakcb) then
      write (this%iout, fmtflow)
    end if
    !
    if (found%iprpak) then
      write (this%iout, '(4x,a)') &
        'LISTS OF '//trim(adjustl(this%text))//' CELLS WILL BE PRINTED.'
    end if
    !
    if (found%iprflow) then
      write (this%iout, '(4x,a)') trim(adjustl(this%text))// &
        ' FLOWS WILL BE PRINTED TO LISTING FILE.'
    end if
    !
    if (found%boundnames) then
      write (this%iout, '(4x,a)') trim(adjustl(this%text))// &
        ' BOUNDARIES HAVE NAMES IN LAST COLUMN.'
    end if
    !
    if (found%auxmultname) then
      write (this%iout, '(4x,a,a)') &
        'AUXILIARY MULTIPLIER NAME: ', sfacauxname
    end if
    !
    if (found%inewton) then
      write (this%iout, '(4x,a)') &
        'NEWTON-RAPHSON method disabled for unconfined cells'
    end if
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' BASE OPTIONS'
  end subroutine log_options

  !> @ brief Source package dimensions from input context
  !<
  subroutine source_dimensions(this)
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy variables
    class(BndExtType), intent(inout) :: this !< BndExtType object
    ! -- local variables
    type(BndExtFoundType) :: found
    !
    ! -- open dimensions logging block
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
      ' BASE DIMENSIONS'
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%maxbound, 'MAXBOUND', this%input_mempath, &
                       found%maxbound)
    !
    write (this%iout, '(4x,a,i7)') 'MAXBOUND = ', this%maxbound
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' BASE DIMENSIONS'
    !
    ! -- verify dimensions were set
    if (this%maxbound <= 0) then
      write (errmsg, '(a)') 'MAXBOUND must be an integer greater than zero.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
  end subroutine source_dimensions

  !> @ brief Update package nodelist
    !!
    !! Convert period updated cellids to node numbers.
    !!
  !<
  subroutine nodelist_update(this)
    ! -- modules
    use SimVariablesModule, only: errmsg
    ! -- dummy
    class(BndExtType) :: this !< BndExtType object
    ! -- local
    integer(I4B), dimension(:), pointer :: cellid
    integer(I4B) :: n, nodeu, noder
    character(len=LINELENGTH) :: nodestr
    !
    ! -- update nodelist
    do n = 1, this%nbound
      !
      ! -- set cellid
      cellid => this%cellid(:, n)
      !
      ! -- ensure cellid is valid, store an error otherwise
      call this%check_cellid(n, cellid, this%dis%mshape, this%dis%ndim)
      !
      ! -- Determine user node number
      if (this%dis%ndim == 1) then
        nodeu = cellid(1)
      elseif (this%dis%ndim == 2) then
        nodeu = get_node(cellid(1), 1, cellid(2), &
                         this%dis%mshape(1), 1, &
                         this%dis%mshape(2))
      else
        nodeu = get_node(cellid(1), cellid(2), cellid(3), &
                         this%dis%mshape(1), &
                         this%dis%mshape(2), &
                         this%dis%mshape(3))
      end if
      !
      ! -- update the nodelist
      if (this%dis%nodes < this%dis%nodesuser) then
        ! -- convert user to reduced node numbers
        noder = this%dis%get_nodenumber(nodeu, 0)
        if (noder <= 0) then
          call this%dis%nodeu_to_string(nodeu, nodestr)
          write (errmsg, *) &
            ' Cell is outside active grid domain: '// &
            trim(adjustl(nodestr))
          call store_error(errmsg)
        end if
        this%nodelist(n) = noder
      else
        this%nodelist(n) = nodeu
      end if
    end do
    !
    ! -- exit if errors were found
    if (count_errors() > 0) then
      write (errmsg, *) count_errors(), ' errors encountered.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
  end subroutine nodelist_update

  !> @ brief Check for valid cellid
  !<
  subroutine check_cellid(this, ii, cellid, mshape, ndim)
    ! -- modules
    use SimVariablesModule, only: errmsg
    ! -- dummy
    class(BndExtType) :: this !< BndExtType object
    ! -- local
    integer(I4B), intent(in) :: ii
    integer(I4B), dimension(:), intent(in) :: cellid !< cellid
    integer(I4B), dimension(:), intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: ndim !< size of mshape
    character(len=20) :: cellstr, mshstr
    character(len=*), parameter :: fmterr = &
      "('List entry ',i0,' contains cellid ',a,' but this cellid is invalid &
      &for model with shape ', a)"
    character(len=*), parameter :: fmtndim1 = &
                                   "('(',i0,')')"
    character(len=*), parameter :: fmtndim2 = &
                                   "('(',i0,',',i0,')')"
    character(len=*), parameter :: fmtndim3 = &
                                   "('(',i0,',',i0,',',i0,')')"
    select case (ndim)
    case (1)
      !
      if (cellid(1) < 1 .or. cellid(1) > mshape(1)) then
        write (cellstr, fmtndim1) cellid(1)
        write (mshstr, fmtndim1) mshape(1)
        write (errmsg, fmterr) ii, trim(adjustl(cellstr)), trim(adjustl(mshstr))
        call store_error(errmsg)
      end if
      !
    case (2)
      !
      if (cellid(1) < 1 .or. cellid(1) > mshape(1) .or. &
          cellid(2) < 1 .or. cellid(2) > mshape(2)) then
        write (cellstr, fmtndim2) cellid(1), cellid(2)
        write (mshstr, fmtndim2) mshape(1), mshape(2)
        write (errmsg, fmterr) ii, trim(adjustl(cellstr)), trim(adjustl(mshstr))
        call store_error(errmsg)
      end if
      !
    case (3)
      !
      if (cellid(1) < 1 .or. cellid(1) > mshape(1) .or. &
          cellid(2) < 1 .or. cellid(2) > mshape(2) .or. &
          cellid(3) < 1 .or. cellid(3) > mshape(3)) then
        write (cellstr, fmtndim3) cellid(1), cellid(2), cellid(3)
        write (mshstr, fmtndim3) mshape(1), mshape(2), mshape(3)
        write (errmsg, fmterr) ii, trim(adjustl(cellstr)), trim(adjustl(mshstr))
        call store_error(errmsg)
      end if
      !
    case default
    end select
  end subroutine check_cellid

  !> @ brief Log package list input
    !!
    !! Log period list based input. This routine requires a package specific
    !! bound_value() routine to report accurate bound values.
    !!
  !<
  subroutine write_list(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENBOUNDNAME, &
                               TABLEFT, TABCENTER, DZERO
    use InputOutputModule, only: ulstlb
    use TableModule, only: TableType, table_cr
    ! -- dummy
    class(BndExtType) :: this !< BndExtType object
    ! -- local
    character(len=10) :: cpos
    character(len=LINELENGTH) :: tag
    character(len=LINELENGTH), allocatable, dimension(:) :: words
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: ipos
    integer(I4B) :: ii, jj, i, j, k, nod
    integer(I4B) :: ldim
    integer(I4B) :: naux
    type(TableType), pointer :: inputtab => null()
    ! -- formats
    character(len=LINELENGTH) :: fmtlstbn
    !
    ! -- Determine sizes
    ldim = this%ncolbnd
    naux = size(this%auxvar, 1)
    !
    ! -- dimension table
    ntabrows = this%nbound
    !
    ! -- start building format statement to parse this%label, which
    !    contains the column headers (except for boundname and auxnames)
    ipos = index(this%listlabel, 'NO.')
    if (ipos /= 0) then
      write (cpos, '(i10)') ipos + 3
      fmtlstbn = '(a'//trim(adjustl(cpos))
    else
      fmtlstbn = '(a7'
    end if
    ! -- sequence number, layer, row, and column.
    if (size(this%dis%mshape) == 3) then
      ntabcols = 4
      fmtlstbn = trim(fmtlstbn)//',a7,a7,a7'
      !
      ! -- sequence number, layer, and cell2d.
    else if (size(this%dis%mshape) == 2) then
      ntabcols = 3
      fmtlstbn = trim(fmtlstbn)//',a7,a7'
      !
      ! -- sequence number and node.
    else
      ntabcols = 2
      fmtlstbn = trim(fmtlstbn)//',a7'
    end if
    !
    ! -- Add fields for non-optional real values
    ntabcols = ntabcols + ldim
    do i = 1, ldim
      fmtlstbn = trim(fmtlstbn)//',a16'
    end do
    !
    ! -- Add field for boundary name
    if (this%inamedbound == 1) then
      ntabcols = ntabcols + 1
      fmtlstbn = trim(fmtlstbn)//',a16'
    end if
    !
    ! -- Add fields for auxiliary variables
    ntabcols = ntabcols + naux
    do i = 1, naux
      fmtlstbn = trim(fmtlstbn)//',a16'
    end do
    fmtlstbn = trim(fmtlstbn)//')'
    !
    ! -- allocate words
    allocate (words(ntabcols))
    !
    ! -- parse this%listlabel into words
    read (this%listlabel, fmtlstbn) (words(i), i=1, ntabcols)
    !
    ! -- initialize the input table object
    call table_cr(inputtab, ' ', ' ')
    call inputtab%table_df(ntabrows, ntabcols, this%iout)
    !
    ! -- add the columns
    ipos = 1
    call inputtab%initialize_column(words(ipos), 10, alignment=TABCENTER)
    !
    ! -- discretization
    do i = 1, size(this%dis%mshape)
      ipos = ipos + 1
      call inputtab%initialize_column(words(ipos), 7, alignment=TABCENTER)
    end do
    !
    ! -- non-optional variables
    do i = 1, ldim
      ipos = ipos + 1
      call inputtab%initialize_column(words(ipos), 16, alignment=TABCENTER)
    end do
    !
    ! -- boundname
    if (this%inamedbound == 1) then
      ipos = ipos + 1
      tag = 'BOUNDNAME'
      call inputtab%initialize_column(tag, LENBOUNDNAME, alignment=TABLEFT)
    end if
    !
    ! -- aux variables
    do i = 1, naux
      call inputtab%initialize_column(this%auxname(i), 16, alignment=TABCENTER)
    end do
    !
    ! -- Write the table
    do ii = 1, this%nbound
      call inputtab%add_term(ii)
      !
      ! -- discretization
      if (size(this%dis%mshape) == 3) then
        nod = this%nodelist(ii)
        call get_ijk(nod, this%dis%mshape(2), this%dis%mshape(3), &
                     this%dis%mshape(1), i, j, k)
        call inputtab%add_term(k)
        call inputtab%add_term(i)
        call inputtab%add_term(j)
      else if (size(this%dis%mshape) == 2) then
        nod = this%nodelist(ii)
        call get_ijk(nod, 1, this%dis%mshape(2), this%dis%mshape(1), i, j, k)
        call inputtab%add_term(k)
        call inputtab%add_term(j)
      else
        nod = this%nodelist(ii)
        call inputtab%add_term(nod)
      end if
      !
      ! -- non-optional variables
      do jj = 1, ldim
        call inputtab%add_term(this%bound_value(jj, ii))
      end do
      !
      ! -- boundname
      if (this%inamedbound == 1) then
        call inputtab%add_term(this%boundname(ii))
      end if
      !
      ! -- aux variables
      do jj = 1, naux
        call inputtab%add_term(this%auxvar(jj, ii))
      end do
    end do
    !
    ! -- deallocate the local variables
    call inputtab%table_da()
    deallocate (inputtab)
    nullify (inputtab)
    deallocate (words)
  end subroutine write_list

  !> @ brief Return a bound value
    !!
    !!  Return a bound value associated with an ncolbnd index
    !!  and row.  This function should be overridden in the
    !!  derived package class.
    !!
  !<
  function bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DNODATA
    ! -- dummy variables
    class(BndExtType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    ! -- override this return value by redefining this
    !    routine in the derived package.
    bndval = DNODATA
  end function bound_value

end module BndExtModule
