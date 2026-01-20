module VirtualGweModelModule
  use KindModule, only: I4B
  use SimStagesModule
  use VirtualBaseModule
  use VirtualDataContainerModule, only: VDC_GWEMODEL_TYPE
  use VirtualModelModule
  use NumericalModelModule, only: NumericalModelType
  implicit none
  private

  public :: add_virtual_gwe_model

  type, extends(VirtualModelType) :: VirtualGweModelType
    ! CND
    type(VirtualIntType), pointer :: cnd_idisp => null()
    type(VirtualDbl1dType), pointer :: cnd_alh => null()
    type(VirtualDbl1dType), pointer :: cnd_alv => null()
    type(VirtualDbl1dType), pointer :: cnd_ath1 => null()
    type(VirtualDbl1dType), pointer :: cnd_ath2 => null()
    type(VirtualDbl1dType), pointer :: cnd_atv => null()
    type(VirtualDbl1dType), pointer :: cnd_ktw => null()
    type(VirtualDbl1dType), pointer :: cnd_kts => null()
    ! FMI
    type(VirtualIntType), pointer :: fmi_igwfspdis => null()
    type(VirtualDbl1dType), pointer :: fmi_gwfhead => null()
    type(VirtualDbl1dType), pointer :: fmi_gwfsat => null()
    type(VirtualDbl2dType), pointer :: fmi_gwfspdis => null()
    type(VirtualDbl1dType), pointer :: fmi_gwfflowja => null()
    ! EST
    type(VirtualDbl1dType), pointer :: est_porosity => null()
    ! GWE Model fields
    type(VirtualIntType), pointer :: incnd => null()
    type(VirtualIntType), pointer :: inest => null()
  contains
    ! public
    procedure :: create => vgwe_create
    procedure :: prepare_stage => vgwe_prepare_stage
    procedure :: destroy => vgwe_destroy
    ! private
    procedure, private :: init_virtual_data
    procedure, private :: allocate_data
    procedure, private :: deallocate_data
  end type VirtualGweModelType

contains

  subroutine add_virtual_gwe_model(model_id, model_name, model)
    use VirtualDataListsModule, only: virtual_model_list
    integer(I4B) :: model_id !< global model id
    character(len=*) :: model_name !< model name
    class(NumericalModelType), pointer :: model !< the actual model (can be null() when remote)
    ! local
    class(VirtualGweModelType), pointer :: virtual_gwe_model
    class(*), pointer :: obj

    allocate (virtual_gwe_model)
    call virtual_gwe_model%create(model_name, model_id, model)

    obj => virtual_gwe_model
    call virtual_model_list%Add(obj)

  end subroutine add_virtual_gwe_model

  subroutine vgwe_create(this, name, id, model)
    class(VirtualGweModelType) :: this
    character(len=*) :: name
    integer(I4B) :: id
    class(NumericalModelType), pointer :: model

    ! create base
    call this%VirtualModelType%create(name, id, model)
    this%container_type = VDC_GWEMODEL_TYPE

    call this%allocate_data()
    call this%init_virtual_data()

  end subroutine vgwe_create

  subroutine init_virtual_data(this)
    class(VirtualGweModelType) :: this

    call this%set(this%cnd_idisp%base(), 'IDISP', 'CND', MAP_ALL_TYPE)
    call this%set(this%cnd_alh%base(), 'ALH', 'CND', MAP_NODE_TYPE)
    call this%set(this%cnd_alv%base(), 'ALV', 'CND', MAP_NODE_TYPE)
    call this%set(this%cnd_ath1%base(), 'ATH1', 'CND', MAP_NODE_TYPE)
    call this%set(this%cnd_ath2%base(), 'ATH2', 'CND', MAP_NODE_TYPE)
    call this%set(this%cnd_atv%base(), 'ATV', 'CND', MAP_NODE_TYPE)
    call this%set(this%cnd_ktw%base(), 'KTW', 'CND', MAP_NODE_TYPE)
    call this%set(this%cnd_kts%base(), 'KTS', 'CND', MAP_NODE_TYPE)
    call this%set(this%fmi_igwfspdis%base(), 'IGWFSPDIS', 'FMI', MAP_ALL_TYPE)
    call this%set(this%fmi_gwfhead%base(), 'GWFHEAD', 'FMI', MAP_NODE_TYPE)
    call this%set(this%fmi_gwfsat%base(), 'GWFSAT', 'FMI', MAP_NODE_TYPE)
    call this%set(this%fmi_gwfspdis%base(), 'GWFSPDIS', 'FMI', MAP_NODE_TYPE)
    call this%set(this%fmi_gwfflowja%base(), 'GWFFLOWJA', 'FMI', MAP_CONN_TYPE)
    call this%set(this%est_porosity%base(), 'POROSITY', 'EST', MAP_NODE_TYPE)
    call this%set(this%incnd%base(), 'INCND', '', MAP_ALL_TYPE)
    call this%set(this%inest%base(), 'INEST', '', MAP_ALL_TYPE)

  end subroutine init_virtual_data

  subroutine vgwe_prepare_stage(this, stage)
    class(VirtualGweModelType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: nr_nodes, nr_conns

    ! prepare base (=numerical) model data items
    call this%VirtualModelType%prepare_stage(stage)

    nr_nodes = 0
    nr_conns = 0

    if (stage == STG_AFT_MDL_DF) then

      call this%map(this%cnd_idisp%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%incnd%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%inest%base(), (/STG_AFT_MDL_DF/))

    else if (stage == STG_BFR_CON_AR) then

      nr_nodes = this%element_maps(MAP_NODE_TYPE)%nr_virt_elems
      nr_conns = this%element_maps(MAP_CONN_TYPE)%nr_virt_elems

      call this%map(this%x%base(), nr_nodes, &
                    (/STG_BFR_CON_AR, STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
      call this%map(this%ibound%base(), nr_nodes, (/STG_BFR_CON_AR/))

      call this%map(this%fmi_igwfspdis%base(), (/STG_BFR_CON_AR/))

      if (this%cnd_idisp%get() > 0) then
        call this%map(this%cnd_alh%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%cnd_alv%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%cnd_ath1%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%cnd_ath2%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%cnd_atv%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%cnd_ktw%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%cnd_kts%base(), nr_nodes, (/STG_BFR_CON_AR/))
      end if

      if (this%incnd%get() > 0 .and. this%inest%get() > 0) then
        call this%map(this%est_porosity%base(), nr_nodes, (/STG_AFT_CON_AR/))
      end if

    else if (stage == STG_AFT_CON_AR) then

      nr_nodes = this%element_maps(MAP_NODE_TYPE)%nr_virt_elems
      nr_conns = this%element_maps(MAP_CONN_TYPE)%nr_virt_elems

      call this%map(this%fmi_gwfhead%base(), nr_nodes, (/STG_BFR_EXG_AD/))
      call this%map(this%fmi_gwfsat%base(), nr_nodes, (/STG_BFR_EXG_AD/))
      if (this%fmi_igwfspdis%get() > 0) then
        call this%map(this%fmi_gwfspdis%base(), 3, nr_nodes, (/STG_BFR_EXG_AD/))
      else
        call this%map(this%fmi_gwfspdis%base(), 3, 0, (/STG_NEVER/))
      end if
      call this%map(this%fmi_gwfflowja%base(), nr_conns, (/STG_BFR_EXG_AD/))

    end if

  end subroutine vgwe_prepare_stage

  subroutine allocate_data(this)
    class(VirtualGweModelType) :: this

    allocate (this%cnd_idisp)
    allocate (this%cnd_alh)
    allocate (this%cnd_alv)
    allocate (this%cnd_ath1)
    allocate (this%cnd_ath2)
    allocate (this%cnd_atv)
    allocate (this%cnd_ktw)
    allocate (this%cnd_kts)
    allocate (this%fmi_igwfspdis)
    allocate (this%fmi_gwfhead)
    allocate (this%fmi_gwfsat)
    allocate (this%fmi_gwfspdis)
    allocate (this%fmi_gwfflowja)
    allocate (this%est_porosity)
    allocate (this%incnd)
    allocate (this%inest)

  end subroutine allocate_data

  subroutine deallocate_data(this)
    class(VirtualGweModelType) :: this

    deallocate (this%cnd_idisp)
    deallocate (this%cnd_alh)
    deallocate (this%cnd_alv)
    deallocate (this%cnd_ath1)
    deallocate (this%cnd_ath2)
    deallocate (this%cnd_atv)
    deallocate (this%cnd_ktw)
    deallocate (this%cnd_kts)
    deallocate (this%fmi_igwfspdis)
    deallocate (this%fmi_gwfhead)
    deallocate (this%fmi_gwfsat)
    deallocate (this%fmi_gwfspdis)
    deallocate (this%fmi_gwfflowja)
    deallocate (this%est_porosity)
    deallocate (this%incnd)
    deallocate (this%inest)

  end subroutine deallocate_data

  subroutine vgwe_destroy(this)
    class(VirtualGweModelType) :: this

    call this%VirtualModelType%destroy()
    call this%deallocate_data()

  end subroutine vgwe_destroy

end module VirtualGweModelModule
