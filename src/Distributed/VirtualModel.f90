module VirtualModelModule
  use VirtualBaseModule
  use VirtualDataContainerModule
  use ConstantsModule, only: LENMEMPATH
  use KindModule, only: I4B, LGP
  use ListModule, only: ListType
  use SimStagesModule
  use NumericalModelModule, only: NumericalModelType
  implicit none
  private

  public :: cast_as_virtual_model
  public :: get_virtual_model_from_list
  public :: get_virtual_model

  interface get_virtual_model
    module procedure get_virtual_model_by_id, &
      get_virtual_model_by_name
  end interface

  type, public, extends(VirtualDataContainerType) :: VirtualModelType
    class(NumericalModelType), pointer :: local_model
    ! CON
    type(VirtualIntType), pointer :: con_ianglex => null()
    type(VirtualInt1dType), pointer :: con_ia => null()
    type(VirtualInt1dType), pointer :: con_ja => null()
    type(VirtualInt1dType), pointer :: con_jas => null()
    type(VirtualInt1dType), pointer :: con_ihc => null()
    type(VirtualDbl1dType), pointer :: con_hwva => null()
    type(VirtualDbl1dType), pointer :: con_cl1 => null()
    type(VirtualDbl1dType), pointer :: con_cl2 => null()
    type(VirtualDbl1dType), pointer :: con_anglex => null()
    ! DIS
    type(VirtualIntType), pointer :: dis_ndim => null()
    type(VirtualIntType), pointer :: dis_nodes => null()
    type(VirtualIntType), pointer :: dis_nodesuser => null()
    type(VirtualInt1dType), pointer :: dis_nodeuser => null()
    type(VirtualIntType), pointer :: dis_nja => null()
    type(VirtualIntType), pointer :: dis_njas => null()
    type(VirtualDblType), pointer :: dis_xorigin => null()
    type(VirtualDblType), pointer :: dis_yorigin => null()
    type(VirtualDblType), pointer :: dis_angrot => null()
    type(VirtualIntType), pointer :: dis_icondir => null()
    type(VirtualDbl1dType), pointer :: dis_xc => null()
    type(VirtualDbl1dType), pointer :: dis_yc => null()
    type(VirtualDbl1dType), pointer :: dis_top => null()
    type(VirtualDbl1dType), pointer :: dis_bot => null()
    type(VirtualDbl1dType), pointer :: dis_area => null()
    ! Numerical Model fields
    type(VirtualIntType), pointer :: moffset => null()
    type(VirtualDbl1dType), pointer :: x => null()
    type(VirtualDbl1dType), pointer :: x_old => null()
    type(VirtualInt1dType), pointer :: ibound => null()
    ! Base Model fields
    type(VirtualIntType), pointer :: idsoln => null()
  contains
    ! public
    procedure :: create => vm_create
    procedure :: prepare_stage => vm_prepare_stage
    procedure :: destroy => vm_destroy
    generic :: operator(==) => eq_virtual_model, eq_numerical_model

    procedure :: dis_get_nodeuser
    procedure :: dis_noder_to_string

    ! private
    procedure, private :: init_virtual_data
    procedure, private :: allocate_data
    procedure, private :: deallocate_data
    procedure, private :: eq_virtual_model
    procedure, private :: eq_numerical_model
  end type VirtualModelType

contains

  subroutine vm_create(this, name, id, model)
    class(VirtualModelType) :: this
    character(len=*) :: name
    integer(I4B) :: id
    class(NumericalModelType), pointer :: model
    ! local
    logical(LGP) :: is_local

    is_local = associated(model)
    call this%VirtualDataContainerType%vdc_create(name, id, is_local)

    this%local_model => model

    call this%allocate_data()
    call this%init_virtual_data()

  end subroutine vm_create

  subroutine init_virtual_data(this)
    class(VirtualModelType) :: this

    ! CON
    call this%set(this%con_ianglex%base(), 'IANGLEX', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_ia%base(), 'IA', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_ja%base(), 'JA', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_jas%base(), 'JAS', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_ihc%base(), 'IHC', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_hwva%base(), 'HWVA', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_cl1%base(), 'CL1', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_cl2%base(), 'CL2', 'CON', MAP_ALL_TYPE)
    call this%set(this%con_anglex%base(), 'ANGLEX', 'CON', MAP_ALL_TYPE)
    ! DIS
    call this%set(this%dis_ndim%base(), 'NDIM', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_nodes%base(), 'NODES', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_nodesuser%base(), 'NODESUSER', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_nodeuser%base(), 'NODEUSER', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_nja%base(), 'NJA', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_njas%base(), 'NJAS', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_icondir%base(), 'ICONDIR', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_xorigin%base(), 'XORIGIN', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_yorigin%base(), 'YORIGIN', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_angrot%base(), 'ANGROT', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_xc%base(), 'XC', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_yc%base(), 'YC', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_top%base(), 'TOP', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_bot%base(), 'BOT', 'DIS', MAP_ALL_TYPE)
    call this%set(this%dis_area%base(), 'AREA', 'DIS', MAP_ALL_TYPE)
    ! Numerical model
    call this%set(this%moffset%base(), 'MOFFSET', '', MAP_ALL_TYPE)
    call this%set(this%x%base(), 'X', '', MAP_NODE_TYPE)
    call this%set(this%x_old%base(), 'XOLD', '', MAP_NODE_TYPE)
    call this%set(this%ibound%base(), 'IBOUND', '', MAP_NODE_TYPE)
    ! Base model
    call this%set(this%idsoln%base(), 'IDSOLN', '', MAP_ALL_TYPE)

  end subroutine init_virtual_data

  subroutine vm_prepare_stage(this, stage)
    class(VirtualModelType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: nodes, nodesuser, nja, njas, icondir
    logical(LGP) :: is_reduced

    if (stage == STG_AFT_MDL_DF) then

      call this%map(this%idsoln%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%con_ianglex%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dis_icondir%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dis_ndim%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dis_nodes%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dis_nodesuser%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dis_nja%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dis_njas%base(), (/STG_AFT_MDL_DF/))

    else if (stage == STG_BFR_EXG_AC) then

      nodes = this%dis_nodes%get()
      nodesuser = this%dis_nodesuser%get()
      is_reduced = (nodes /= nodesuser)
      call this%map(this%moffset%base(), (/STG_BFR_EXG_AC/))
      if (is_reduced) then
        call this%map(this%dis_nodeuser%base(), nodes, (/STG_BFR_EXG_AC/))
      else
        ! no reduction, zero sized array, never synchronize
        call this%map(this%dis_nodeuser%base(), 0, (/STG_NEVER/))
      end if

    else if (stage == STG_BFR_CON_DF) then

      nodes = this%dis_nodes%get()
      nja = this%dis_nja%get()
      njas = this%dis_njas%get()
      icondir = this%dis_icondir%get()
      ! DIS
      call this%map(this%dis_xorigin%base(), (/STG_BFR_CON_DF/))
      call this%map(this%dis_yorigin%base(), (/STG_BFR_CON_DF/))
      call this%map(this%dis_angrot%base(), (/STG_BFR_CON_DF/))
      if (icondir > 0) then
        call this%map(this%dis_xc%base(), nodes, (/STG_BFR_CON_DF/))
        call this%map(this%dis_yc%base(), nodes, (/STG_BFR_CON_DF/))
      else
        call this%map(this%dis_xc%base(), 0, (/STG_NEVER/))
        call this%map(this%dis_yc%base(), 0, (/STG_NEVER/))
      end if
      call this%map(this%dis_top%base(), nodes, (/STG_BFR_CON_DF/))
      call this%map(this%dis_bot%base(), nodes, (/STG_BFR_CON_DF/))
      call this%map(this%dis_area%base(), nodes, (/STG_BFR_CON_DF/))
      ! CON
      call this%map(this%con_ia%base(), nodes + 1, (/STG_BFR_CON_DF/))
      call this%map(this%con_ja%base(), nja, (/STG_BFR_CON_DF/))
      call this%map(this%con_jas%base(), nja, (/STG_BFR_CON_DF/))
      call this%map(this%con_ihc%base(), njas, (/STG_BFR_CON_DF/))
      call this%map(this%con_hwva%base(), njas, (/STG_BFR_CON_DF/))
      call this%map(this%con_cl1%base(), njas, (/STG_BFR_CON_DF/))
      call this%map(this%con_cl2%base(), njas, (/STG_BFR_CON_DF/))
      if (this%con_ianglex%get() > 0) then
        call this%map(this%con_anglex%base(), njas, (/STG_BFR_CON_DF/))
      else
        call this%map(this%con_anglex%base(), 0, (/STG_NEVER/))
      end if

    end if

  end subroutine vm_prepare_stage

  !> @brief Get user node number from reduced number
  !<
  function dis_get_nodeuser(this, node_reduced) result(node_user)
    class(VirtualModelType) :: this !< this virtual model
    integer(I4B), intent(in) :: node_reduced !< the reduced node number
    integer(I4B) :: node_user !< the returned user node number

    if (this%dis_nodes%get() < this%dis_nodesuser%get()) then
      node_user = this%dis_nodeuser%get(node_reduced)
    else
      node_user = node_reduced
    end if

  end function dis_get_nodeuser

  subroutine dis_noder_to_string(this, node_reduced, node_str)
    class(VirtualModelType) :: this !< this virtual model
    integer(I4B), intent(in) :: node_reduced !< reduced node number
    character(len=*), intent(inout) :: node_str !< the string representative of the user node number
    ! local
    character(len=11) :: nr_str

    if (this%is_local) then
      call this%local_model%dis%noder_to_string(node_reduced, node_str)
    else
      ! for now this will look like: (102r)
      write (nr_str, '(i0)') node_reduced
      node_str = '('//trim(adjustl(nr_str))//'r)'
    end if

  end subroutine dis_noder_to_string

  subroutine vm_destroy(this)
    class(VirtualModelType) :: this

    call this%VirtualDataContainerType%destroy()
    call this%deallocate_data()

  end subroutine vm_destroy

  subroutine allocate_data(this)
    class(VirtualModelType) :: this

    allocate (this%con_ianglex)
    allocate (this%con_ia)
    allocate (this%con_ja)
    allocate (this%con_jas)
    allocate (this%con_ihc)
    allocate (this%con_hwva)
    allocate (this%con_cl1)
    allocate (this%con_cl2)
    allocate (this%con_anglex)
    allocate (this%dis_ndim)
    allocate (this%dis_nodes)
    allocate (this%dis_nodesuser)
    allocate (this%dis_nodeuser)
    allocate (this%dis_nja)
    allocate (this%dis_njas)
    allocate (this%dis_xorigin)
    allocate (this%dis_yorigin)
    allocate (this%dis_angrot)
    allocate (this%dis_icondir)
    allocate (this%dis_xc)
    allocate (this%dis_yc)
    allocate (this%dis_top)
    allocate (this%dis_bot)
    allocate (this%dis_area)
    allocate (this%moffset)
    allocate (this%x)
    allocate (this%x_old)
    allocate (this%ibound)
    allocate (this%idsoln)

  end subroutine allocate_data

  subroutine deallocate_data(this)
    class(VirtualModelType) :: this

    ! CON
    deallocate (this%con_ianglex)
    deallocate (this%con_ia)
    deallocate (this%con_ja)
    deallocate (this%con_jas)
    deallocate (this%con_ihc)
    deallocate (this%con_hwva)
    deallocate (this%con_cl1)
    deallocate (this%con_cl2)
    deallocate (this%con_anglex)
    ! DIS
    deallocate (this%dis_ndim)
    deallocate (this%dis_nodes)
    deallocate (this%dis_nodesuser)
    deallocate (this%dis_nodeuser)
    deallocate (this%dis_nja)
    deallocate (this%dis_njas)
    deallocate (this%dis_xorigin)
    deallocate (this%dis_yorigin)
    deallocate (this%dis_angrot)
    deallocate (this%dis_icondir)
    deallocate (this%dis_xc)
    deallocate (this%dis_yc)
    deallocate (this%dis_top)
    deallocate (this%dis_bot)
    deallocate (this%dis_area)
    ! Numerical model
    deallocate (this%moffset)
    deallocate (this%x)
    deallocate (this%x_old)
    deallocate (this%ibound)
    ! Base model
    deallocate (this%idsoln)

  end subroutine deallocate_data

  function get_virtual_model_from_list(model_list, idx) result(v_model)
    type(ListType) :: model_list
    integer(I4B) :: idx
    class(VirtualModelType), pointer :: v_model
    ! local
    class(*), pointer :: obj_ptr

    obj_ptr => model_list%GetItem(idx)
    v_model => cast_as_virtual_model(obj_ptr)
  end function get_virtual_model_from_list

  function cast_as_virtual_model(obj_ptr) result(v_model)
    class(*), pointer :: obj_ptr
    class(VirtualModelType), pointer :: v_model

    v_model => null()
    select type (obj_ptr)
    class is (VirtualModelType)
      v_model => obj_ptr
    end select

  end function cast_as_virtual_model

  function eq_virtual_model(this, v_model) result(is_equal)
    class(VirtualModelType), intent(in) :: this
    class(VirtualModelType), intent(in) :: v_model
    logical(LGP) :: is_equal

    is_equal = (this%id == v_model%id)

  end function eq_virtual_model

  function eq_numerical_model(this, num_model) result(is_equal)
    class(VirtualModelType), intent(in) :: this
    class(NumericalModelType), intent(in) :: num_model
    logical(LGP) :: is_equal

    is_equal = (this%id == num_model%id)

  end function eq_numerical_model

  !> @brief Returns a virtual model with the specified id
  !< from the global list, or null
  function get_virtual_model_by_id(model_id) result(virtual_model)
    use VirtualDataListsModule, only: virtual_model_list
    integer(I4B) :: model_id
    class(VirtualModelType), pointer :: virtual_model
    ! local
    integer(I4B) :: i
    class(*), pointer :: vm

    virtual_model => null()
    do i = 1, virtual_model_list%Count()
      vm => virtual_model_list%GetItem(i)
      select type (vm)
      class is (VirtualModelType)
        if (vm%id == model_id) then
          virtual_model => vm
          return
        end if
      end select
    end do

  end function get_virtual_model_by_id

  !> @brief Returns a virtual model with the specified name
  !< from the global list, or null
  function get_virtual_model_by_name(model_name) result(virtual_model)
    use VirtualDataListsModule, only: virtual_model_list
    character(len=*) :: model_name
    class(VirtualModelType), pointer :: virtual_model
    ! local
    integer(I4B) :: i
    class(*), pointer :: vm

    virtual_model => null()
    do i = 1, virtual_model_list%Count()
      vm => virtual_model_list%GetItem(i)
      select type (vm)
      class is (VirtualModelType)
        if (vm%name == model_name) then
          virtual_model => vm
          return
        end if
      end select
    end do

  end function get_virtual_model_by_name

end module VirtualModelModule
