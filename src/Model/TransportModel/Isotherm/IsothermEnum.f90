module IsothermEnumModule
  use KindModule, only: I4B

  implicit none

  integer(I4B), parameter :: SORPTION_OFF = 0 !< Sorption is inactive (default)
  integer(I4B), parameter :: SORPTION_LINEAR = 1 !< Linear sorption between aqueous and solid phases
  integer(I4B), parameter :: SORPTION_FREUND = 2 !< Freundlich sorption between aqueous and solid phases
  integer(I4B), parameter :: SORPTION_LANG = 3 !< Langmuir sorption between aqueous and solid phases

end module IsothermEnumModule
