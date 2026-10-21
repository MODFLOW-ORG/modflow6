!> @brief Source and check the static UZF cell properties
!<
submodule(UzfModule) UzfModuleInput
contains

  !> @brief Source PACKAGEDATA block from input context
  !<
  module procedure uzf_source_packagedata
  ! -- modules
  use SimModule, only: store_error, store_error_filename, count_errors
  ! -- local
  integer(I4B), dimension(:), pointer, contiguous :: ifno => null()
  integer(I4B), dimension(:, :), pointer, contiguous :: cellid_idm => null()
  integer(I4B), dimension(:), pointer, contiguous :: landflag => null()
  integer(I4B), dimension(:), pointer, contiguous :: ivertcon => null()
  real(DP), dimension(:), pointer, contiguous :: surfdep => null()
  real(DP), dimension(:), pointer, contiguous :: vks => null()
  real(DP), dimension(:), pointer, contiguous :: thtr => null()
  real(DP), dimension(:), pointer, contiguous :: thts => null()
  real(DP), dimension(:), pointer, contiguous :: thti => null()
  real(DP), dimension(:), pointer, contiguous :: eps => null()
  type(CharacterStringType), dimension(:), pointer, contiguous :: &
    bndname => null()
  character(len=LINELENGTH) :: cellid_str
  integer(I4B) :: ierr
  integer(I4B) :: i
  integer(I4B) :: n !< uzf cell number
  integer(I4B) :: node !< gwf node number
  integer(I4B) :: j, jcol
  integer(I4B), dimension(:), allocatable :: rowmaxnnz
  type(sparsematrix) :: sparse
  integer(I4B), dimension(:), allocatable :: nboundchk
  !
  write (this%iout, '(/1x,3a)') 'PROCESSING ', trim(adjustl(this%text)), &
    ' PACKAGEDATA'
  !
  ! -- get PACKAGEDATA arrays
  call mem_setptr(ifno, 'PACKAGEDATA_IFNO', this%input_mempath)
  call mem_setptr(cellid_idm, 'CELLID', this%input_mempath)
  call mem_setptr(landflag, 'LANDFLAG', this%input_mempath)
  call mem_setptr(ivertcon, 'IVERTCON', this%input_mempath)
  call mem_setptr(surfdep, 'SURFDEP', this%input_mempath)
  call mem_setptr(vks, 'VKS', this%input_mempath)
  call mem_setptr(thtr, 'THTR', this%input_mempath)
  call mem_setptr(thts, 'THTS', this%input_mempath)
  call mem_setptr(thti, 'THTI', this%input_mempath)
  call mem_setptr(eps, 'EPS', this%input_mempath)
  if (this%inamedbound /= 0) then
    call mem_setptr(bndname, 'BOUNDNAME', this%input_mempath)
  end if
  !
  ! -- allocate local check arrays
  allocate (rowmaxnnz(this%dis%nodes))
  do node = 1, this%dis%nodes
    rowmaxnnz(node) = 0
  end do
  allocate (nboundchk(this%nodes))
  do n = 1, this%nodes
    nboundchk(n) = 0
  end do
  !
  ! -- loop over packagedata entries
  do i = 1, this%nodes
    n = this%validate_ifno(ifno(i), this%nodes, nboundchk, &
                           'IUZNO', 'PACKAGEDATA')
    if (n == 0) cycle
    !
    ! -- convert CELLID to node number
    write (cellid_str, '(10(i0,1x))') cellid_idm(1:this%dis%ndim, i)
    node = this%dis%noder_from_cellid(cellid_str, this%inunit, this%iout)
    this%igwfnode(n) = node
    rowmaxnnz(node) = rowmaxnnz(node) + 1
    !
    ! -- validate landflag
    if (landflag(i) < 0 .OR. landflag(i) > 1) then
      write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
        'LANDFLAG for uzf cell', n, &
        'must be 0 or 1 (specified value is', landflag(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- validate ivertcon
    if (ivertcon(i) < 0 .OR. ivertcon(i) > this%nodes) then
      write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
        'IVERTCON for uzf cell', n, &
        'must be 0 or less than NUZFCELLS (specified value is', &
        ivertcon(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- validate surfdep
    if (surfdep(i) <= DZERO .and. landflag(i) > 0) then
      write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
        'SURFDEP for uzf cell', n, &
        'must be greater than 0 (specified value is', surfdep(i), ').'
      call store_error(errmsg)
    end if
    if (surfdep(i) >= this%dis%top(node) - this%dis%bot(node)) then
      write (errmsg, '(a,1x,i0,1x,a)') &
        'SURFDEP for uzf cell', n, &
        'cannot be greater than the cell thickness.'
      call store_error(errmsg)
    end if
    !
    ! -- validate vks
    if (vks(i) <= DZERO) then
      write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
        'VKS for uzf cell', n, &
        'must be greater than 0 (specified value ia', vks(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- validate thtr
    if (thtr(i) <= DZERO) then
      write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
        'THTR for uzf cell', n, &
        'must be greater than 0 (specified value is', thtr(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- validate thts
    if (thts(i) <= thtr(i)) then
      write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
        'THTS for uzf cell', n, &
        'must be greater than THTR (specified value is', thts(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- validate thti
    if (thti(i) < thtr(i) .OR. thti(i) > thts(i)) then
      write (errmsg, '(a,1x,i0,1x,a,1x,a,1x,g0,a)') &
        'THTI for uzf cell', n, &
        'must be greater than or equal to THTR AND less than THTS', &
        '(specified value is', thti(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- validate eps
    if (eps(i) < 3.5 .OR. eps(i) > 14) then
      write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
        'EPSILON for uzf cell', n, &
        'must be between 3.5 and 14.0 (specified value is', eps(i), ').'
      call store_error(errmsg)
    end if
    !
    ! -- boundname
    if (this%inamedbound == 1) then
      this%uzfname(n) = bndname(i)
    end if
    !
    ! -- set data if no errors
    if (count_errors() == 0) then
      call this%uzfobj%setdata(n, this%dis%area(node), &
                               this%dis%top(node), this%dis%bot(node), &
                               surfdep(i), vks(i), thtr(i), thts(i), &
                               thti(i), eps(i), this%ntrail_input, &
                               landflag(i), ivertcon(i))
      if (ivertcon(i) > 0) then
        this%iuzf2uzf = 1
      end if
    end if
  end do
  write (this%iout, '(1x,3a)') &
    'END OF ', trim(adjustl(this%text)), ' PACKAGEDATA'
  !
  ! -- check for duplicate or missing uzf cells
  call this%report_ifno_coverage(nboundchk, this%nodes, 'uzf cell', &
                                 'PACKAGEDATA')
  if (count_errors() > 0) then
    call store_error_filename(this%input_fname)
  end if
  !
  ! -- setup sparse for connectivity used to identify multiple uzf cells
  !    per GWF model cell
  call sparse%init(this%dis%nodes, this%dis%nodes, rowmaxnnz)
  do n = 1, this%nodes
    node = this%igwfnode(n)
    call sparse%addconnection(node, n, 1)
  end do
  !
  ! -- create ia and ja from sparse
  call sparse%filliaja(this%ia, this%ja, ierr)
  !
  ! -- set imaxcellcnt
  do i = 1, this%dis%nodes
    jcol = 0
    do j = this%ia(i), this%ia(i + 1) - 1
      jcol = jcol + 1
    end do
    if (jcol > this%imaxcellcnt) then
      this%imaxcellcnt = jcol
    end if
  end do
  !
  if (this%imaxcellcnt > 1 .and. this%iauxmultcol < 1) then
    call this%check_cell_area()
  end if
  !
  deallocate (rowmaxnnz)
  deallocate (nboundchk)
  !
  ! -- release input context packagedata memory
  call memorystore_release('PACKAGEDATA_IFNO', this%input_mempath)
  call memorystore_release('CELLID', this%input_mempath)
  call memorystore_release('LANDFLAG', this%input_mempath)
  call memorystore_release('IVERTCON', this%input_mempath)
  call memorystore_release('SURFDEP', this%input_mempath)
  call memorystore_release('VKS', this%input_mempath)
  call memorystore_release('THTR', this%input_mempath)
  call memorystore_release('THTS', this%input_mempath)
  call memorystore_release('THTI', this%input_mempath)
  call memorystore_release('EPS', this%input_mempath)
  call memorystore_release('BOUNDNAME', this%input_mempath)
  end procedure uzf_source_packagedata

  !> @brief Read UZF cell properties and set them for UZFCellGroup type
  !<
  module procedure print_cell_properties
  ! -- local
  character(len=20) :: cellid
  character(len=LINELENGTH) :: title
  character(len=LINELENGTH) :: tag
  integer(I4B) :: ntabrows
  integer(I4B) :: ntabcols
  integer(I4B) :: n
  integer(I4B) :: node
  !
  ntabrows = this%nodes
  ntabcols = 10
  if (this%inamedbound == 1) then
    ntabcols = ntabcols + 1
  end if
  !
  ! -- initialize table and define columns
  title = trim(adjustl(this%text))//' PACKAGE ('// &
          trim(adjustl(this%packName))//') STATIC UZF CELL DATA'
  call table_cr(this%inputtab, this%packName, title)
  call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
  tag = 'NUMBER'
  call this%inputtab%initialize_column(tag, 10)
  tag = 'CELLID'
  call this%inputtab%initialize_column(tag, 20, alignment=TABLEFT)
  tag = 'LANDFLAG'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'IVERTCON'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'SURFDEP'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'VKS'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'THTR'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'THTS'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'THTI'
  call this%inputtab%initialize_column(tag, 12)
  tag = 'EPS'
  call this%inputtab%initialize_column(tag, 12)
  if (this%inamedbound == 1) then
    tag = 'BOUNDNAME'
    call this%inputtab%initialize_column(tag, LENBOUNDNAME, alignment=TABLEFT)
  end if
  !
  ! -- write data for each cell
  do n = 1, this%nodes
    !
    ! -- get cellid
    node = this%igwfnode(n)
    if (node > 0) then
      call this%dis%noder_to_string(node, cellid)
    else
      cellid = 'none'
    end if
    !
    ! -- add data
    call this%inputtab%add_term(n)
    call this%inputtab%add_term(cellid)
    call this%inputtab%add_term(this%uzfobj%landflag(n))
    call this%inputtab%add_term(this%uzfobj%cell_below(n))
    call this%inputtab%add_term(this%uzfobj%surfdep(n))
    call this%inputtab%add_term(this%uzfobj%vks(n))
    call this%inputtab%add_term(this%uzfobj%theta_res(n))
    call this%inputtab%add_term(this%uzfobj%theta_sat(n))
    call this%inputtab%add_term(this%uzfobj%theta_init(n))
    call this%inputtab%add_term(this%uzfobj%bc_eps(n))
    if (this%inamedbound == 1) then
      call this%inputtab%add_term(this%uzfname(n))
    end if
  end do
  end procedure print_cell_properties

  !> @brief Check UZF cell areas
  !<
  module procedure check_cell_area
  ! -- modules
  use InputOutputModule, only: urword
  use SimModule, only: store_error, count_errors, store_error_filename
  ! -- local
  character(len=16) :: cuzf
  character(len=20) :: cellid
  character(len=LINELENGTH) :: cuzfcells
  integer(I4B) :: n
  integer(I4B) :: n2
  integer(I4B) :: j
  integer(I4B) :: node
  integer(I4B) :: i0
  integer(I4B) :: i1
  real(DP) :: area
  real(DP) :: area2
  real(DP) :: sumarea
  real(DP) :: cellarea
  real(DP) :: d
  !
  do n = 1, this%nodes
    !
    ! -- Initialize variables
    n2 = this%uzfobj%cell_below(n)
    area = this%uzfobj%uzfarea(n)
    !
    ! Create pointer to object below
    if (n2 > 0) then
      area2 = this%uzfobj%uzfarea(n2)
      d = abs(area - area2)
      if (d > DEM6) then
        write (errmsg, '(2(a,1x,g0,1x,a,1x,i0,1x),a)') &
          'UZF cell area (', area, ') for cell ', n, &
          'does not equal uzf cell area (', area2, ') for cell ', n2, '.'
        call store_error(errmsg)
      end if
    end if
  end do
  !
  ! -- check that the area of uzf cells in a GWF cell is less than or equal
  !    to the GWF cell area
  do node = 1, this%dis%nodes
    i0 = this%ia(node)
    i1 = this%ia(node + 1)
    ! -- skip gwf cells with no UZF cells
    if ((i1 - i0) < 1) cycle
    sumarea = DZERO
    cellarea = DZERO
    cuzfcells = ''
    do j = i0, i1 - 1
      n = this%ja(j)
      write (cuzf, '(i0)') n
      cuzfcells = trim(adjustl(cuzfcells))//' '//trim(adjustl(cuzf))
      sumarea = sumarea + this%uzfobj%uzfarea(n)
      cellarea = this%uzfobj%cellarea(n)
    end do
    ! -- calculate the difference between the sum of UZF areas and GWF cell area
    d = sumarea - cellarea
    if (d > DEM6) then
      call this%dis%noder_to_string(node, cellid)
      write (errmsg, '(a,1x,g0,1x,a,1x,g0,1x,a,1x,a,1x,a,a,a)') &
        'Total uzf cell area (', sumarea, &
        ') exceeds the gwf cell area (', cellarea, ') of cell', cellid, &
        'which includes uzf cell(s): ', trim(adjustl(cuzfcells)), '.'
      call store_error(errmsg)
    end if
  end do
  !
  ! -- terminate if errors were encountered
  if (count_errors() > 0) then
    call store_error_filename(this%input_fname)
  end if
  end procedure check_cell_area

end submodule UzfModuleInput
