module VirtualGwtModelModule
  use KindModule, only: I4B
  use SimStagesModule
  use VirtualBaseModule
  use VirtualDataContainerModule, only: VDC_GWTMODEL_TYPE
  use VirtualModelModule
  use NumericalModelModule, only: NumericalModelType
  implicit none
  private

  public :: add_virtual_gwt_model

  type, extends(VirtualModelType) :: VirtualGwtModelType
    ! DSP
    type(VirtualIntType), pointer :: dsp_idiffc => null()
    type(VirtualIntType), pointer :: dsp_idisp => null()
    type(VirtualDbl1dType), pointer :: dsp_diffc => null()
    type(VirtualDbl1dType), pointer :: dsp_alh => null()
    type(VirtualDbl1dType), pointer :: dsp_alv => null()
    type(VirtualDbl1dType), pointer :: dsp_ath1 => null()
    type(VirtualDbl1dType), pointer :: dsp_ath2 => null()
    type(VirtualDbl1dType), pointer :: dsp_atv => null()
    ! FMI
    type(VirtualIntType), pointer :: fmi_igwfspdis => null()
    type(VirtualDbl1dType), pointer :: fmi_gwfhead => null()
    type(VirtualDbl1dType), pointer :: fmi_gwfsat => null()
    type(VirtualDbl2dType), pointer :: fmi_gwfspdis => null()
    type(VirtualDbl1dType), pointer :: fmi_gwfflowja => null()
    ! MST
    type(VirtualDbl1dType), pointer :: mst_thetam => null()
    ! GWT Model fields
    type(VirtualIntType), pointer :: indsp => null()
    type(VirtualIntType), pointer :: inmst => null()
  contains
    ! public
    procedure :: create => vgwt_create
    procedure :: prepare_stage => vgwt_prepare_stage
    procedure :: destroy => vgwt_destroy
    ! private
    procedure, private :: init_virtual_data
    procedure, private :: allocate_data
    procedure, private :: deallocate_data
  end type VirtualGwtModelType

contains

  subroutine add_virtual_gwt_model(model_id, model_name, model)
    use VirtualDataListsModule, only: virtual_model_list
    integer(I4B) :: model_id !< global model id
    character(len=*) :: model_name !< model name
    class(NumericalModelType), pointer :: model !< the actual model (can be null() when remote)
    ! local
    class(VirtualGwtModelType), pointer :: virtual_gwt_model
    class(*), pointer :: obj

    allocate (virtual_gwt_model)
    call virtual_gwt_model%create(model_name, model_id, model)

    obj => virtual_gwt_model
    call virtual_model_list%Add(obj)

  end subroutine add_virtual_gwt_model

  subroutine vgwt_create(this, name, id, model)
    class(VirtualGwtModelType) :: this
    character(len=*) :: name
    integer(I4B) :: id
    class(NumericalModelType), pointer :: model

    ! create base
    call this%VirtualModelType%create(name, id, model)
    this%container_type = VDC_GWTMODEL_TYPE

    call this%allocate_data()
    call this%init_virtual_data()

  end subroutine vgwt_create

  subroutine init_virtual_data(this)
    class(VirtualGwtModelType) :: this

    call this%set(this%dsp_idiffc%base(), 'IDIFFC', 'DSP', MAP_ALL_TYPE)
    call this%set(this%dsp_idisp%base(), 'IDISP', 'DSP', MAP_ALL_TYPE)
    call this%set(this%dsp_diffc%base(), 'DIFFC', 'DSP', MAP_NODE_TYPE)
    call this%set(this%dsp_alh%base(), 'ALH', 'DSP', MAP_NODE_TYPE)
    call this%set(this%dsp_alv%base(), 'ALV', 'DSP', MAP_NODE_TYPE)
    call this%set(this%dsp_ath1%base(), 'ATH1', 'DSP', MAP_NODE_TYPE)
    call this%set(this%dsp_ath2%base(), 'ATH2', 'DSP', MAP_NODE_TYPE)
    call this%set(this%dsp_atv%base(), 'ATV', 'DSP', MAP_NODE_TYPE)
    call this%set(this%fmi_igwfspdis%base(), 'IGWFSPDIS', 'FMI', MAP_ALL_TYPE)
    call this%set(this%fmi_gwfhead%base(), 'GWFHEAD', 'FMI', MAP_NODE_TYPE)
    call this%set(this%fmi_gwfsat%base(), 'GWFSAT', 'FMI', MAP_NODE_TYPE)
    call this%set(this%fmi_gwfspdis%base(), 'GWFSPDIS', 'FMI', MAP_NODE_TYPE)
    call this%set(this%fmi_gwfflowja%base(), 'GWFFLOWJA', 'FMI', MAP_CONN_TYPE)
    call this%set(this%mst_thetam%base(), 'THETAM', 'MST', MAP_NODE_TYPE)
    call this%set(this%indsp%base(), 'INDSP', '', MAP_ALL_TYPE)
    call this%set(this%inmst%base(), 'INMST', '', MAP_ALL_TYPE)

  end subroutine init_virtual_data

  subroutine vgwt_prepare_stage(this, stage)
    class(VirtualGwtModelType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: nr_nodes, nr_conns

    ! prepare base (=numerical) model data items
    call this%VirtualModelType%prepare_stage(stage)

    nr_nodes = 0
    nr_conns = 0

    if (stage == STG_AFT_MDL_DF) then

      call this%map(this%dsp_idiffc%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%dsp_idisp%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%indsp%base(), (/STG_AFT_MDL_DF/))
      call this%map(this%inmst%base(), (/STG_AFT_MDL_DF/))

    else if (stage == STG_BFR_CON_AR) then

      nr_nodes = this%element_maps(MAP_NODE_TYPE)%nr_virt_elems
      nr_conns = this%element_maps(MAP_CONN_TYPE)%nr_virt_elems

      call this%map(this%fmi_igwfspdis%base(), (/STG_BFR_CON_AR/))
      call this%map(this%x%base(), nr_nodes, &
                    (/STG_BFR_CON_AR, STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
      call this%map(this%ibound%base(), nr_nodes, (/STG_BFR_CON_AR/))

      if (this%dsp_idiffc%get() > 0) then
        call this%map(this%dsp_diffc%base(), nr_nodes, (/STG_BFR_CON_AR/))
      end if

      if (this%dsp_idisp%get() > 0) then
        call this%map(this%dsp_alh%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%dsp_alv%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%dsp_ath1%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%dsp_ath2%base(), nr_nodes, (/STG_BFR_CON_AR/))
        call this%map(this%dsp_atv%base(), nr_nodes, (/STG_BFR_CON_AR/))
      end if

      if (this%indsp%get() > 0 .and. this%inmst%get() > 0) then
        call this%map(this%mst_thetam%base(), nr_nodes, (/STG_AFT_CON_AR/))
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

  end subroutine vgwt_prepare_stage

  subroutine allocate_data(this)
    class(VirtualGwtModelType) :: this

    allocate (this%dsp_idiffc)
    allocate (this%dsp_idisp)
    allocate (this%dsp_diffc)
    allocate (this%dsp_alh)
    allocate (this%dsp_alv)
    allocate (this%dsp_ath1)
    allocate (this%dsp_ath2)
    allocate (this%dsp_atv)
    allocate (this%fmi_igwfspdis)
    allocate (this%fmi_gwfhead)
    allocate (this%fmi_gwfsat)
    allocate (this%fmi_gwfspdis)
    allocate (this%fmi_gwfflowja)
    allocate (this%mst_thetam)
    allocate (this%indsp)
    allocate (this%inmst)

  end subroutine allocate_data

  subroutine deallocate_data(this)
    class(VirtualGwtModelType) :: this

    deallocate (this%dsp_idiffc)
    deallocate (this%dsp_idisp)
    deallocate (this%dsp_diffc)
    deallocate (this%dsp_alh)
    deallocate (this%dsp_alv)
    deallocate (this%dsp_ath1)
    deallocate (this%dsp_ath2)
    deallocate (this%dsp_atv)
    deallocate (this%fmi_igwfspdis)
    deallocate (this%fmi_gwfhead)
    deallocate (this%fmi_gwfsat)
    deallocate (this%fmi_gwfspdis)
    deallocate (this%fmi_gwfflowja)
    deallocate (this%mst_thetam)
    deallocate (this%indsp)
    deallocate (this%inmst)

  end subroutine deallocate_data

  subroutine vgwt_destroy(this)
    class(VirtualGwtModelType) :: this

    call this%VirtualModelType%destroy()
    call this%deallocate_data()

  end subroutine vgwt_destroy

end module VirtualGwtModelModule
