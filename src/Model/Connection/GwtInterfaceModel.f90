module GwtInterfaceModelModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DONE
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_reallocate
  use MemoryHelperModule, only: create_mem_path
  use NumericalModelModule, only: NumericalModelType
  use GwtModule, only: GwtModelType, CastAsGwtModel
  use DisuModule, only: disu_cr, CastAsDisuType
  use TspFmiModule, only: fmi_cr, TspFmiType
  use TspAdvModule, only: adv_cr, TspAdvType
  use TspAdvOptionsModule, only: TspAdvOptionsType
  use GwtDspModule, only: dsp_cr, GwtDspType
  use GwtDspOptionsModule, only: GwtDspOptionsType
  use GwtMstModule, only: mst_cr
  use TspObsModule, only: tsp_obs_cr
  use GridConnectionModule

  implicit none
  private

  !> The GWT Interface Model is a utility to calculate the solution's exchange
  !! coefficients from the interface between a GWT model and its GWT neighbors.
  !! The interface model itself will not be part of the solution, it is not
  !! being solved.
  !<
  type, public, extends(GwtModelType) :: GwtInterfaceModelType

    integer(i4B), pointer :: iAdvScheme => null() !< the advection scheme: 0 = up, 1 = central, 2 = tvd
    integer(i4B), pointer :: ixt3d => null() !< xt3d setting: 0 = off, 1 = lhs, 2 = rhs
    real(DP), pointer :: ieqnsclfac => null() !< governing eqn scaling factor: 1: GWT, >1: GWE

    class(GridConnectionType), pointer :: gridConnection => null() !< The grid connection class will provide the interface grid
    class(GwtModelType), private, pointer :: owner => null() !< the real GWT model for which the exchange coefficients
                                                             !! are calculated with this interface model

  contains
    procedure, pass(this) :: gwtifmod_cr
    procedure :: model_df => gwtifmod_df
    procedure :: model_ar => gwtifmod_ar
    procedure :: model_da => gwtifmod_da
    procedure, public :: allocate_fmi
    procedure :: allocate_scalars
  end type GwtInterfaceModelType

contains

  !> @brief Create the interface model, analogously to what
  !< happens in gwt_cr
  subroutine gwtifmod_cr(this, name, iout, gridConn)
    ! -- dummy
    class(GwtInterfaceModelType) :: this !< the GWT interface model
    character(len=*), intent(in) :: name !< the interface model's name
    integer(I4B), intent(in) :: iout !< the output unit
    class(GridConnectionType), pointer, intent(in) :: gridConn !< the grid connection data for creating a DISU
    ! local
    class(*), pointer :: modelPtr
    integer(I4B), target :: inobs
    integer(I4B) :: adv_unit, dsp_unit
    !
    this%memoryPath = create_mem_path(name)
    call this%allocate_scalars(name)
    !
    ! defaults
    this%iAdvScheme = 0
    this%ixt3d = 0
    this%ieqnsclfac = DONE
    !
    this%iout = iout
    this%gridConnection => gridConn
    modelPtr => gridConn%model
    this%owner => CastAsGwtModel(modelPtr)
    !
    inobs = 0
    adv_unit = 0
    dsp_unit = 0
    if (this%owner%inadv > 0) then
      this%inadv = huge(1_I4B)
      adv_unit = huge(1_I4B)
    end if
    if (this%owner%indsp > 0) then
      this%indsp = huge(1_I4B)
      dsp_unit = huge(1_I4B)
    end if
    !
    ! create dis and packages
    call disu_cr(this%dis, this%name, '', -1, this%iout)
    call fmi_cr(this%fmi, this%name, 0, this%iout, this%ieqnsclfac, &
                this%depvartype)
    call adv_cr(this%adv, this%name, '', adv_unit, this%iout, this%fmi, &
                this%ieqnsclfac)
    call dsp_cr(this%dsp, this%name, '', -dsp_unit, this%iout, this%fmi)
    call tsp_obs_cr(this%obs, inobs, this%depvartype)
  end subroutine gwtifmod_cr

  !> @brief Allocate scalars associated with the interface model object
  !<
  subroutine allocate_scalars(this, modelname)
    ! -- dummy
    class(GwtInterfaceModelType) :: this !< the GWT interface model
    character(len=*), intent(in) :: modelname !< the model name
    !
    call this%GwtModelType%allocate_scalars(modelname)
    !
    call mem_allocate(this%iAdvScheme, 'ADVSCHEME', this%memoryPath)
    call mem_allocate(this%ixt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%ieqnsclfac, 'IEQNSCLFAC', this%memoryPath)
  end subroutine allocate_scalars

  !> @brief Allocate a Flow Model Interface (FMI) object for the interface model
  !<
  subroutine allocate_fmi(this)
    ! -- dummy
    class(GwtInterfaceModelType) :: this !< the GWT interface model
    !
    call mem_allocate(this%fmi%gwfflowja, this%nja, 'GWFFLOWJA', &
                      this%fmi%memoryPath)
    call mem_allocate(this%fmi%gwfhead, this%neq, 'GWFHEAD', &
                      this%fmi%memoryPath)
    call mem_allocate(this%fmi%gwfsat, this%neq, 'GWFSAT', &
                      this%fmi%memoryPath)
    call mem_allocate(this%fmi%gwfspdis, 3, this%neq, 'GWFSPDIS', &
                      this%fmi%memoryPath)
  end subroutine allocate_fmi

  !> @brief Define the GWT interface model
  !<
  subroutine gwtifmod_df(this)
    ! -- dummy
    class(GwtInterfaceModelType) :: this !< the GWT interface model
    ! -- local
    class(*), pointer :: disPtr
    type(TspAdvOptionsType) :: adv_options
    type(GwtDspOptionsType) :: dsp_options
    !
    this%moffset = 0
    adv_options%iAdvScheme = this%iAdvScheme
    dsp_options%ixt3d = this%ixt3d
    !
    ! define DISU
    disPtr => this%dis
    call this%gridConnection%getDiscretization(CastAsDisuType(disPtr))
    call this%fmi%fmi_df(this%dis, 1)
    !
    if (this%inadv > 0) then
      call this%adv%adv_df(adv_options)
    end if
    if (this%indsp > 0) then
      this%dsp%idiffc = this%owner%dsp%idiffc
      this%dsp%idisp = this%owner%dsp%idisp
      call this%dsp%dsp_df(this%dis, dsp_options)
      if (this%dsp%idiffc > 0) then
        call mem_reallocate(this%dsp%diffc, this%dis%nodes, 'DIFFC', &
                            trim(this%dsp%memoryPath))
      end if
      if (this%dsp%idisp > 0) then
        call mem_reallocate(this%dsp%alh, this%dis%nodes, 'ALH', &
                            trim(this%dsp%memoryPath))
        call mem_reallocate(this%dsp%alv, this%dis%nodes, 'ALV', &
                            trim(this%dsp%memoryPath))
        call mem_reallocate(this%dsp%ath1, this%dis%nodes, 'ATH1', &
                            trim(this%dsp%memoryPath))
        call mem_reallocate(this%dsp%ath2, this%dis%nodes, 'ATH2', &
                            trim(this%dsp%memoryPath))
        call mem_reallocate(this%dsp%atv, this%dis%nodes, 'ATV', &
                            trim(this%dsp%memoryPath))
      end if
      allocate (this%mst)
      call mem_allocate(this%mst%thetam, this%dis%nodes, &
                        'THETAM', create_mem_path(this%name, 'MST'))
    end if
    !
    ! assign or point model members to dis members
    this%neq = this%dis%nodes
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja
    !
    ! allocate model arrays, now that neq and nja are assigned
    call this%allocate_arrays()
  end subroutine gwtifmod_df

  !> @brief Override allocate and read the GWT interface model and its packages
  !! so that we can create stuff from memory instead of input files
  !<
  subroutine gwtifmod_ar(this)
    ! -- dummy
    class(GwtInterfaceModelType) :: this !< the GWT interface model
    !
    call this%fmi%fmi_ar(this%ibound)
    if (this%inadv > 0) then
      call this%adv%adv_ar(this%dis, this%ibound)
    end if
    if (this%indsp > 0) then
      call this%dsp%dsp_ar(this%ibound, this%mst%thetam)
    end if
  end subroutine gwtifmod_ar

  !> @brief Clean up resources
  !<
  subroutine gwtifmod_da(this)
    ! -- dummy
    class(GwtInterfaceModelType) :: this !< the GWT interface model
    !

    ! this
    call mem_deallocate(this%iAdvScheme)
    call mem_deallocate(this%ixt3d)
    call mem_deallocate(this%ieqnsclfac)
    !
    ! gwt packages
    call this%dis%dis_da()
    call this%fmi%fmi_da()
    call this%adv%adv_da()
    call this%dsp%dsp_da()
    !
    deallocate (this%dis)
    deallocate (this%fmi)
    deallocate (this%adv)
    deallocate (this%dsp)
    !
    if (associated(this%mst)) then
      call mem_deallocate(this%mst%thetam)
      deallocate (this%mst)
    end if
    !
    ! gwt scalars
    call mem_deallocate(this%inic)
    call mem_deallocate(this%infmi)
    call mem_deallocate(this%inadv)
    call mem_deallocate(this%indsp)
    call mem_deallocate(this%inssm)
    call mem_deallocate(this%inmst)
    call mem_deallocate(this%inmvt)
    call mem_deallocate(this%inoc)
    call mem_deallocate(this%inobs)
    call mem_deallocate(this%eqnsclfac)
    !
    ! base
    call this%NumericalModelType%model_da()
  end subroutine gwtifmod_da

end module GwtInterfaceModelModule
