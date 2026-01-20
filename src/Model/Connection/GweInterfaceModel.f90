module GweInterfaceModelModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DONE
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_reallocate
  use MemoryHelperModule, only: create_mem_path
  use NumericalModelModule, only: NumericalModelType
  use GweModule, only: GweModelType, CastAsGweModel
  use DisuModule, only: disu_cr, CastAsDisuType
  use TspFmiModule, only: fmi_cr, TspFmiType
  use TspAdvModule, only: adv_cr, TspAdvType
  use TspAdvOptionsModule, only: TspAdvOptionsType
  use GweCndModule, only: cnd_cr, GweCndType
  use GweCndOptionsModule, only: GweCndOptionsType
  use GweEstModule, only: est_cr
  use TspObsModule, only: tsp_obs_cr
  use GridConnectionModule
  use GweInputDataModule, only: GweInputDataType

  implicit none
  private

  !> The GWE Interface Model is a utility to calculate the solution's exchange
  !! coefficients from the interface between a GWE model and its GWE neighbors.
  !! The interface model itself will not be part of the solution, it is not
  !! being solved.
  !<
  type, public, extends(GweModelType) :: GweInterfaceModelType

    integer(i4B), pointer :: iAdvScheme => null() !< the advection scheme: 0 = up, 1 = central, 2 = tvd
    integer(i4B), pointer :: ixt3d => null() !< xt3d setting: 0 = off, 1 = lhs, 2 = rhs
    real(DP), pointer :: ieqnsclfac => null() !< governing eqn scaling factor: 1: GWT, >1: GWE

    class(GridConnectionType), pointer :: gridConnection => null() !< The grid connection class will provide the interface grid
    class(GweModelType), private, pointer :: owner => null() !< the real GWE model for which the exchange coefficients
                                                             !! are calculated with this interface model

    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< to be filled with EST porosity

  contains

    procedure, pass(this) :: gweifmod_cr
    procedure :: model_df => gweifmod_df
    procedure :: model_ar => gweifmod_ar
    procedure :: model_da => gweifmod_da
    procedure, public :: allocate_fmi
    procedure :: allocate_scalars
  end type GweInterfaceModelType

contains

  !> @brief Create the interface model, analogously to what
  !< happens in gwe_cr
  subroutine gweifmod_cr(this, name, iout, gridConn)
    ! -- modules
    use GweInputDataModule, only: gweshared_dat_cr
    ! -- dummy
    class(GweInterfaceModelType) :: this !< the GWE interface model
    character(len=*), intent(in) :: name !< the interface model's name
    integer(I4B), intent(in) :: iout !< the output unit
    class(GridConnectionType), pointer, intent(in) :: gridConn !< the grid connection data for creating a DISU
    ! -- local
    class(*), pointer :: modelPtr
    integer(I4B), target :: inobs
    integer(I4B) :: adv_unit, cnd_unit
    !
    this%memoryPath = create_mem_path(name)
    call this%allocate_scalars(name)
    !
    ! -- Instantiate shared data container
    call gweshared_dat_cr(this%gwecommon)
    !
    ! -- Defaults
    this%iAdvScheme = 0
    this%ixt3d = 0
    this%ieqnsclfac = DONE
    !
    this%iout = iout
    this%gridConnection => gridConn
    modelPtr => gridConn%model
    this%owner => CastAsGweModel(modelPtr)
    !
    inobs = 0
    adv_unit = 0
    cnd_unit = 0
    if (this%owner%inadv > 0) then
      this%inadv = huge(1_I4B)
      adv_unit = huge(1_I4B)
    end if
    if (this%owner%incnd > 0) then
      this%incnd = huge(1_I4B)
      cnd_unit = huge(1_I4B)
    end if
    !
    ! -- Create dis and packages
    call disu_cr(this%dis, this%name, '', -1, this%iout)
    call fmi_cr(this%fmi, this%name, '', 0, this%iout, this%ieqnsclfac, &
                this%depvartype)
    call adv_cr(this%adv, this%name, '', adv_unit, this%iout, this%fmi, &
                this%ieqnsclfac)
    call cnd_cr(this%cnd, this%name, '', -cnd_unit, this%iout, this%fmi, &
                this%ieqnsclfac, this%gwecommon)
    call tsp_obs_cr(this%obs, inobs, this%depvartype)
  end subroutine gweifmod_cr

  !> @brief Allocate scalars associated with the interface model object
  !<
  subroutine allocate_scalars(this, modelname)
    ! -- dummy
    class(GweInterfaceModelType) :: this !< the GWE interface model
    character(len=*), intent(in) :: modelname !< the model name
    !
    call this%GweModelType%allocate_scalars(modelname)
    !
    call mem_allocate(this%iAdvScheme, 'ADVSCHEME', this%memoryPath)
    call mem_allocate(this%ixt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%ieqnsclfac, 'IEQNSCLFAC', this%memoryPath)
  end subroutine allocate_scalars

  !> @brief Allocate a Flow Model Interface (FMI) object for the interface model
  !<
  subroutine allocate_fmi(this)
    ! -- dummy
    class(GweInterfaceModelType) :: this !< the GWT interface model
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

  !> @brief Define the GWE interface model
  !<
  subroutine gweifmod_df(this)
    ! -- dummy
    class(GweInterfaceModelType) :: this !< the GWE interface model
    ! -- local
    class(*), pointer :: disPtr
    type(TspAdvOptionsType) :: adv_options
    type(GweCndOptionsType) :: cnd_options
    !
    this%moffset = 0
    adv_options%iAdvScheme = this%iAdvScheme
    cnd_options%ixt3d = this%ixt3d
    !
    ! -- Define DISU
    disPtr => this%dis
    call this%gridConnection%getDiscretization(CastAsDisuType(disPtr))
    call this%fmi%fmi_df(this%dis, 0)
    !
    if (this%inadv > 0) then
      call this%adv%adv_df(adv_options)
    end if
    !
    if (this%incnd > 0) then
      this%cnd%idisp = this%owner%cnd%idisp
      call this%cnd%cnd_df(this%dis, cnd_options)
      !
      if (this%cnd%idisp > 0) then
        call mem_reallocate(this%cnd%alh, this%dis%nodes, 'ALH', &
                            trim(this%cnd%memoryPath))
        call mem_reallocate(this%cnd%alv, this%dis%nodes, 'ALV', &
                            trim(this%cnd%memoryPath))
        call mem_reallocate(this%cnd%ath1, this%dis%nodes, 'ATH1', &
                            trim(this%cnd%memoryPath))
        call mem_reallocate(this%cnd%ath2, this%dis%nodes, 'ATH2', &
                            trim(this%cnd%memoryPath))
        call mem_reallocate(this%cnd%atv, this%dis%nodes, 'ATV', &
                            trim(this%cnd%memoryPath))
        call mem_reallocate(this%cnd%ktw, this%dis%nodes, 'KTW', &
                            trim(this%cnd%memoryPath))
        call mem_reallocate(this%cnd%kts, this%dis%nodes, 'KTS', &
                            trim(this%cnd%memoryPath))
      end if
      allocate (this%est)
      call mem_allocate(this%est%porosity, this%dis%nodes, &
                        'POROSITY', create_mem_path(this%name, 'EST'))
    end if
    !
    ! -- Assign or point model members to dis members
    this%neq = this%dis%nodes
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja
    !
    ! -- Allocate model arrays, now that neq and nja are assigned
    call this%allocate_arrays()
  end subroutine gweifmod_df

  !> @brief Override allocate and read the GWE interface model and its packages
  !! so that we can create stuff from memory instead of input files
  !<
  subroutine gweifmod_ar(this)
    ! -- dummy
    class(GweInterfaceModelType) :: this !< the GWE interface model
    !
    call this%fmi%fmi_ar(this%ibound)
    if (this%inadv > 0) then
      call this%adv%adv_ar(this%dis, this%ibound)
    end if
    if (this%incnd > 0) then
      call this%cnd%cnd_ar(this%ibound, this%est%porosity)
    end if
  end subroutine gweifmod_ar

  !> @brief Clean up resources
  !<
  subroutine gweifmod_da(this)
    ! -- dummy
    class(GweInterfaceModelType) :: this !< the GWE interface model
    !
    ! -- members specified in the interface model input for use by the
    !    interface model
    call mem_deallocate(this%iAdvScheme)
    call mem_deallocate(this%ixt3d)
    call mem_deallocate(this%ieqnsclfac)
    !
    ! -- gwe packages
    call this%dis%dis_da()
    call this%fmi%fmi_da()
    call this%adv%adv_da()
    call this%cnd%cnd_da()
    !
    deallocate (this%dis)
    deallocate (this%fmi)
    deallocate (this%adv)
    deallocate (this%cnd)
    !
    if (associated(this%est)) then
      call mem_deallocate(this%est%porosity)
      deallocate (this%est)
    end if
    !
    ! -- gwe scalars
    call mem_deallocate(this%inic)
    call mem_deallocate(this%infmi)
    call mem_deallocate(this%inadv)
    call mem_deallocate(this%incnd)
    call mem_deallocate(this%inssm)
    call mem_deallocate(this%inest)
    call mem_deallocate(this%inmvt)
    call mem_deallocate(this%inoc)
    call mem_deallocate(this%inobs)
    call mem_deallocate(this%eqnsclfac)
    !
    ! -- base
    call this%NumericalModelType%model_da()
  end subroutine gweifmod_da

end module GweInterfaceModelModule
