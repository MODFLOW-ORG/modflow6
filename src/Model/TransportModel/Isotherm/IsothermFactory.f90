module IsothermFactoryModule

  use KindModule, only: I4B, DP
  use SimModule, only: store_error

  use IsothermEnumModule
  use IsothermInterfaceModule, only: IsothermType
  use LinearIsothermsModule, only: LinearIsothermType
  use FreundlichIsothermModule, only: FreundlichIsothermType
  use LangmuirIsothermModule, only: LangmuirIsothermType

  implicit none
  private
  public :: create_isotherm

contains

  !> @brief Create an isotherm object based on type and parameters.
  !>
  !> Returns a pointer to a concrete `IsothermType` (or `null()` if sorption
  !> is off). Uses `isotherm_type` to select Linear, Freundlich, or Langmuir,
  !> passing `distcoef` and `sp2` as required by the chosen model.
  !<
  function create_isotherm(isotherm_type, distcoef, sp2) result(isotherm)
    ! -- result
    class(IsothermType), pointer :: isotherm !< allocated concrete isotherm or null()
    ! -- dummy
    integer(I4B), intent(in) :: isotherm_type !< enumerator from `IsothermEnumModule`
    real(DP), dimension(:), pointer, contiguous :: distcoef !< primary coefficient (Kd, Kf, or Kl)
    real(DP), dimension(:), pointer, contiguous :: sp2 !< secondary parameter (a for Freundlich, Sbar for Langmuir)

    select case (isotherm_type)
    case (SORPTION_OFF)
      nullify (isotherm)
    case (SORPTION_LINEAR)
      allocate (isotherm, source=LinearIsothermType(distcoef))
    case (SORPTION_FREUND)
      allocate (isotherm, source=FreundlichIsothermType(distcoef, sp2))
    case (SORPTION_LANG)
      allocate (isotherm, source=LangmuirIsothermType(distcoef, sp2))
    case default
      call store_error('Sorption type not implemented.')
    end select

  end function create_isotherm
end module IsothermFactoryModule
