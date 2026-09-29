!> @brief Read and check the static UZF cell properties
!<
submodule(UzfModule) UzfModuleInput
contains

  !> @brief Read UZF cell properties and set them for UzfCellGroup type
  !<
  module procedure read_cell_properties
  ! -- modules
  use InputOutputModule, only: urword
  use SimModule, only: store_error, count_errors
  ! -- local
  character(len=LINELENGTH) :: cellid
  integer(I4B) :: ierr
  integer(I4B) :: n !< uzf cell number
  integer(I4B) :: j
  integer(I4B) :: node !< gwf node number
  integer(I4B) :: jcol
  !
  logical :: isfound, endOfBlock
  integer(I4B) :: landflag
  integer(I4B) :: ivertcon
  real(DP) :: surfdep, vks, thtr, thts, thti, eps, hgwf
  integer(I4B), dimension(:), allocatable :: rowmaxnnz
  type(sparsematrix) :: sparse
  integer(I4B), dimension(:), allocatable :: nboundchk
  !
  ! -- allocate space for node counter and initialize
  allocate (rowmaxnnz(this%dis%nodes))
  do node = 1, this%dis%nodes
    rowmaxnnz(node) = 0
  end do
  !
  ! -- allocate space for local variables
  allocate (nboundchk(this%nodes))
  do n = 1, this%nodes
    nboundchk(n) = 0
  end do
  !
  ! -- initialize variables
  landflag = 0
  ivertcon = 0
  surfdep = DZERO
  vks = DZERO
  thtr = DZERO
  thts = DZERO
  thti = DZERO
  eps = DZERO
  hgwf = DZERO
  !
  ! -- get uzf properties block
  call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                            supportOpenClose=.true.)
  !
  ! -- parse locations block if detected
  if (isfound) then
    write (this%iout, '(/1x,3a)') 'PROCESSING ', trim(adjustl(this%text)), &
      ' PACKAGEDATA'
    do
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      !
      ! -- get uzf cell number
      n = this%parser%GetInteger()

      if (n < 1 .or. n > this%nodes) then
        write (errmsg, '(2(a,1x),i0,a)') &
          'IUZNO must be greater than 0 and less than', &
          'or equal to', this%nodes, '.'
        call store_error(errmsg)
        cycle
      end if
      !
      ! -- increment nboundchk
      nboundchk(n) = nboundchk(n) + 1
      !
      ! -- store the reduced gwf nodenumber in igwfnode
      call this%parser%GetCellid(this%dis%ndim, cellid)
      node = this%dis%noder_from_cellid(cellid, &
                                        this%parser%iuactive, this%iout)
      this%igwfnode(n) = node
      rowmaxnnz(node) = rowmaxnnz(node) + 1
      !
      ! -- landflag
      landflag = this%parser%GetInteger()
      if (landflag < 0 .OR. landflag > 1) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
          'LANDFLAG for uzf cell', n, &
          'must be 0 or 1 (specified value is', landflag, ').'
        call store_error(errmsg)
      end if
      !
      ! -- ivertcon
      ivertcon = this%parser%GetInteger()
      if (ivertcon < 0 .OR. ivertcon > this%nodes) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
          'IVERTCON for uzf cell', n, &
          'must be 0 or less than NUZFCELLS (specified value is', &
          ivertcon, ').'
        call store_error(errmsg)
      end if
      !
      ! -- surfdep
      surfdep = this%parser%GetDouble()
      if (surfdep <= DZERO .and. landflag > 0) then !need to check for cell thickness
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'SURFDEP for uzf cell', n, &
          'must be greater than 0 (specified value is', surfdep, ').'
        call store_error(errmsg)
      end if
      if (surfdep >= this%dis%top(node) - this%dis%bot(node)) then
        write (errmsg, '(a,1x,i0,1x,a)') &
          'SURFDEP for uzf cell', n, &
          'cannot be greater than the cell thickness.'
        call store_error(errmsg)
      end if
      !
      ! -- vks
      vks = this%parser%GetDouble()
      if (vks <= DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'VKS for uzf cell', n, &
          'must be greater than 0 (specified value ia', vks, ').'
        call store_error(errmsg)
      end if
      !
      ! -- thtr
      thtr = this%parser%GetDouble()
      if (thtr <= DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'THTR for uzf cell', n, &
          'must be greater than 0 (specified value is', thtr, ').'
        call store_error(errmsg)
      end if
      !
      ! -- thts
      thts = this%parser%GetDouble()
      if (thts <= thtr) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'THTS for uzf cell', n, &
          'must be greater than THTR (specified value is', thts, ').'
        call store_error(errmsg)
      end if
      !
      ! -- thti
      thti = this%parser%GetDouble()
      if (thti < thtr .OR. thti > thts) then
        write (errmsg, '(a,1x,i0,1x,a,1x,a,1x,g0,a)') &
          'THTI for uzf cell', n, &
          'must be greater than or equal to THTR AND less than THTS', &
          '(specified value is', thti, ').'
        call store_error(errmsg)
      end if
      !
      ! -- eps
      eps = this%parser%GetDouble()
      if (eps < 3.5 .OR. eps > 14) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'EPSILON for uzf cell', n, &
          'must be between 3.5 and 14.0 (specified value is', eps, ').'
        call store_error(errmsg)
      end if
      !
      ! -- boundname
      if (this%inamedbound == 1) then
        call this%parser%GetStringCaps(this%uzfname(n))
      end if
      !
      ! -- set data if there are no data errors
      if (count_errors() == 0) then
        node = this%igwfnode(n)
        call this%uzfobj%setdata(n, this%dis%area(node), this%dis%top(node), &
                                 this%dis%bot(node), surfdep, vks, thtr, thts, &
                                 thti, eps, this%ntrail_input, landflag, &
                                 ivertcon)
        if (ivertcon > 0) then
          this%iuzf2uzf = 1
        end if
      end if
      !
    end do
    write (this%iout, '(1x,3a)') &
      'END OF ', trim(adjustl(this%text)), ' PACKAGEDATA'
  else
    call store_error('Required packagedata block not found.')
  end if
  !
  ! -- check for duplicate or missing uzf cells
  do n = 1, this%nodes
    if (nboundchk(n) == 0) then
      write (errmsg, '(a,1x,i0,a)') &
        'No data specified for uzf cell', n, '.'
      call store_error(errmsg)
    else if (nboundchk(n) > 1) then
      write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
        'Data for uzf cell', n, 'specified', nboundchk(n), 'times.'
      call store_error(errmsg)
    end if
  end do
  !
  ! -- write summary of UZF cell property error messages
  if (count_errors() > 0) then
    call this%parser%StoreErrorUnit()
  end if
  !
  ! -- setup sparse for connectivity used to identify multiple uzf cells per
  !    GWF model cell
  call sparse%init(this%dis%nodes, this%dis%nodes, rowmaxnnz)
  ! --
  do n = 1, this%nodes
    node = this%igwfnode(n)
    call sparse%addconnection(node, n, 1)
  end do
  !
  ! -- create ia and ja from sparse
  call sparse%filliaja(this%ia, this%ja, ierr)
  !
  ! -- set imaxcellcnt
  do node = 1, this%dis%nodes
    jcol = 0
    do j = this%ia(node), this%ia(node + 1) - 1
      jcol = jcol + 1
    end do
    if (jcol > this%imaxcellcnt) then
      this%imaxcellcnt = jcol
    end if
  end do
  !
  ! -- do an initial evaluation of the sum of uzfarea relative to the
  !    GWF cell area in the case that there is more than one UZF object
  !    in a GWF cell and a auxmult value is not being applied to the
  !    calculate the UZF cell area from the GWF cell area.
  if (this%imaxcellcnt > 1 .and. this%iauxmultcol < 1) then
    call this%check_cell_area()
  end if
  !
  ! -- deallocate local variables
  deallocate (rowmaxnnz)
  deallocate (nboundchk)
  end procedure read_cell_properties

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
  use SimModule, only: store_error, count_errors
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
    call this%parser%StoreErrorUnit()
  end if
  end procedure check_cell_area

end submodule UzfModuleInput
