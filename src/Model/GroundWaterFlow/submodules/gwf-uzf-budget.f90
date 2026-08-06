!> @brief Build and fill the UZF budget object
!<
submodule(UzfModule) UzfModuleBudget
contains

  !> @brief Set up the budget object that stores all the uzf flows
  !!
  !! The terms listed here must correspond in number and order to the ones
  !! listed in the uzf_fill_budobj routine
  !<
  module procedure uzf_setup_budobj
  ! -- modules
  use ConstantsModule, only: LENBUDTXT
  ! -- local
  integer(I4B) :: nbudterm
  integer(I4B) :: maxlist, naux
  integer(I4B) :: idx
  integer(I4B) :: nlen
  integer(I4B) :: n, n1, n2 !< uzf cell numbers
  integer(I4B) :: node !< gwf node number
  integer(I4B) :: ivertflag
  real(DP) :: q
  character(len=LENBUDTXT) :: text
  character(len=LENBUDTXT), dimension(1) :: auxtxt
  !
  nlen = 0
  do n = 1, this%nodes
    ivertflag = this%uzfobj%cell_below(n)
    if (ivertflag > 0) then
      nlen = nlen + 1
    end if
  end do
  !
  ! -- Determine the number of uzf budget terms. These are fixed for
  !    the simulation and cannot change.  This includes FLOW-JA-FACE
  !    so they can be written to the binary budget files, but these internal
  !    flows are not included as part of the budget table.
  nbudterm = 4
  if (nlen > 0) nbudterm = nbudterm + 1
  if (this%ietflag /= 0) nbudterm = nbudterm + 1
  if (this%imover == 1) nbudterm = nbudterm + 2
  if (this%naux > 0) nbudterm = nbudterm + 1
  !
  ! -- set up budobj
  call budgetobject_cr(this%budobj, this%packName)
  call this%budobj%budgetobject_df(this%maxbound, nbudterm, 0, 0, &
                                   ibudcsv=this%ibudcsv)
  idx = 0
  !
  ! -- Go through and set up each budget term
  text = '    FLOW-JA-FACE'
  if (nlen > 0) then
    idx = idx + 1
    maxlist = nlen * 2
    naux = 1
    auxtxt(1) = '       FLOW-AREA'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux, auxtxt, ordered_id1=.false.)
    !
    ! -- store connectivity
    call this%budobj%budterm(idx)%reset(nlen * 2)
    q = DZERO
    do n = 1, this%nodes
      ivertflag = this%uzfobj%cell_below(n)
      if (ivertflag > 0) then
        n1 = n
        n2 = ivertflag
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n2, n1, -q)
      end if
    end do
  end if
  !
  ! --
  text = '             GWF'
  idx = idx + 1
  maxlist = this%nodes
  naux = 1
  auxtxt(1) = '       FLOW-AREA'
  call this%budobj%budterm(idx)%initialize(text, &
                                           this%name_model, &
                                           this%packName, &
                                           this%name_model, &
                                           this%name_model, &
                                           maxlist, .false., .true., &
                                           naux, auxtxt)
  call this%budobj%budterm(idx)%reset(this%nodes)
  q = DZERO
  do n = 1, this%nodes
    node = this%igwfnode(n)
    this%qauxcbc(1) = this%uzfobj%uzfarea(n)
    call this%budobj%budterm(idx)%update_term(n, node, q, this%qauxcbc)
  end do
  !
  ! --
  text = '    INFILTRATION'
  idx = idx + 1
  maxlist = this%nodes
  naux = 0
  call this%budobj%budterm(idx)%initialize(text, &
                                           this%name_model, &
                                           this%packName, &
                                           this%name_model, &
                                           this%packName, &
                                           maxlist, .false., .false., &
                                           naux)
  !
  ! --
  text = '         REJ-INF'
  idx = idx + 1
  maxlist = this%nodes
  naux = 0
  call this%budobj%budterm(idx)%initialize(text, &
                                           this%name_model, &
                                           this%packName, &
                                           this%name_model, &
                                           this%packName, &
                                           maxlist, .false., .false., &
                                           naux)
  !
  ! --
  text = '            UZET'
  if (this%ietflag /= 0) then
    idx = idx + 1
    maxlist = this%maxbound
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
  end if
  !
  ! --
  text = '         STORAGE'
  idx = idx + 1
  maxlist = this%nodes
  naux = 1
  auxtxt(1) = '          VOLUME'
  call this%budobj%budterm(idx)%initialize(text, &
                                           this%name_model, &
                                           this%packName, &
                                           this%name_model, &
                                           this%packName, &
                                           maxlist, .false., .false., &
                                           naux, auxtxt)
  !
  ! --
  if (this%imover == 1) then
    !
    ! --
    text = '        FROM-MVR'
    idx = idx + 1
    maxlist = this%nodes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '  REJ-INF-TO-MVR'
    idx = idx + 1
    maxlist = this%nodes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
  end if
  !
  ! --
  naux = this%naux
  if (naux > 0) then
    !
    ! --
    text = '       AUXILIARY'
    idx = idx + 1
    maxlist = this%maxbound
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux, this%auxname)
  end if
  !
  ! -- if uzf flow for each reach are written to the listing file
  if (this%iprflow /= 0) then
    call this%budobj%flowtable_df(this%iout, cellids='GWF')
  end if
  end procedure uzf_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  module procedure uzf_fill_budobj
  ! -- local
  integer(I4B) :: naux
  integer(I4B) :: nlen
  integer(I4B) :: ivertflag
  integer(I4B) :: n, n1, n2 !< uzf cell numbers
  integer(I4B) :: node !< gwf node number
  integer(I4B) :: idx
  real(DP) :: q
  real(DP) :: a
  real(DP) :: top
  real(DP) :: bot
  real(DP) :: thick
  real(DP) :: fm
  real(DP) :: v
  !
  idx = 0
  !
  ! -- FLOW JA FACE
  nlen = 0
  do n = 1, this%nodes
    ivertflag = this%uzfobj%cell_below(n)
    if (ivertflag > 0) then
      nlen = nlen + 1
    end if
  end do
  if (nlen > 0) then
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(nlen * 2)
    do n = 1, this%nodes
      ivertflag = this%uzfobj%cell_below(n)
      if (ivertflag > 0) then
        a = this%uzfobj%uzfarea(n)
        q = this%uzfobj%surf_infil_below(n) * a
        this%qauxcbc(1) = a
        if (q > DZERO) then
          q = -q
        end if
        n1 = n
        n2 = ivertflag
        call this%budobj%budterm(idx)%update_term(n1, n2, q, this%qauxcbc)
        call this%budobj%budterm(idx)%update_term(n2, n1, -q, this%qauxcbc)
      end if
    end do
  end if
  !
  ! -- GWF (LEAKAGE)
  idx = idx + 1
  call this%budobj%budterm(idx)%reset(this%nodes)
  do n = 1, this%nodes
    this%qauxcbc(1) = this%uzfobj%uzfarea(n)
    node = this%igwfnode(n)
    q = -this%rch(n)
    call this%budobj%budterm(idx)%update_term(n, node, q, this%qauxcbc)
  end do
  !
  ! -- INFILTRATION
  idx = idx + 1
  call this%budobj%budterm(idx)%reset(this%nodes)
  do n = 1, this%nodes
    q = this%appliedinf(n)
    call this%budobj%budterm(idx)%update_term(n, n, q)
  end do
  !
  ! -- REJECTED INFILTRATION
  idx = idx + 1
  call this%budobj%budterm(idx)%reset(this%nodes)
  do n = 1, this%nodes
    q = this%rejinf(n)
    if (q > DZERO) then
      q = -q
    end if
    call this%budobj%budterm(idx)%update_term(n, n, q)
  end do
  !
  ! -- UNSATURATED EVT
  if (this%ietflag /= 0) then
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = this%uzet(n)
      if (q > DZERO) then
        q = -q
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
  end if
  !
  ! -- STORAGE
  idx = idx + 1
  call this%budobj%budterm(idx)%reset(this%nodes)
  do n = 1, this%nodes
    q = -this%qsto(n)
    top = this%uzfobj%celtop(n)
    bot = this%uzfobj%water_table(n)
    thick = top - bot
    if (thick > DZERO) then
      fm = thick * (this%wcnew(n) - this%uzfobj%theta_res(n))
      v = fm * this%uzfobj%uzfarea(n)
    else
      v = DZERO
    end if
    ! -- save mobile water volume into aux variable
    this%qauxcbc(1) = v
    call this%budobj%budterm(idx)%update_term(n, n, q, this%qauxcbc)
  end do
  !
  ! -- MOVER
  if (this%imover == 1) then
    !
    ! -- FROM MOVER
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = this%pakmvrobj%get_qfrommvr(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- REJ-INF-TO-MVR
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = this%rejinftomvr(n)
      if (q > DZERO) then
        q = -q
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do

  end if
  !
  ! -- AUXILIARY VARIABLES
  naux = this%naux
  if (naux > 0) then
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = DZERO
      call this%budobj%budterm(idx)%update_term(n, n, q, this%auxvar(:, n))
    end do
  end if
  !
  ! --Terms are filled, now accumulate them for this time step
  call this%budobj%accumulate_terms()
  end procedure uzf_fill_budobj

end submodule UzfModuleBudget
