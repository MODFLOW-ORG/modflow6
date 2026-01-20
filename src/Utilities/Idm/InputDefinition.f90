!> @brief This module contains the InputDefinitionModule
!!
!! This module contains type definitions that represent
!! descriptions of input from modflow 6 definition files.
!!
!<
module InputDefinitionModule

  use KindModule, only: LGP
  use ConstantsModule, only: LINELENGTH, LENVARNAME, LENCOMPONENTNAME

  implicit none
  private
  public :: InputParamDefinitionType, &
            InputBlockDefinitionType

  !> @brief Input parameter definition type
  !!
  !! This type is used to store information for
  !! each modflow input record
  !!
  !<
  type InputParamDefinitionType
    character(len=LENCOMPONENTNAME) :: component_type = '' !< type of component, e.g. GWF
    character(len=LENCOMPONENTNAME) :: subcomponent_type = '' !< type of subcomponent, e.g. CHD
    character(len=LINELENGTH) :: blockname = '' !< input block, e.g. DiMENSIONS
    character(len=LINELENGTH) :: tagname = '' !< parameter user tag name
    character(len=LENVARNAME) :: mf6varname = '' !< parameter internal managed memory name
    character(len=LINELENGTH) :: datatype = '' !< parameter data type
    character(len=LINELENGTH) :: shape = '' !< shape of data type
    character(len=LINELENGTH) :: longname = '' !< description of variable
    logical(LGP) :: required = .false. !< is the parameter required
    logical(LGP) :: developmode = .false. !< is the parameter developmode
    logical(LGP) :: in_record = .false. !< is the parameter within an input record
    logical(LGP) :: preserve_case = .false. !< should string case be preserved
    logical(LGP) :: layered = .false. !< does the parameter support a layered read
    logical(LGP) :: timeseries = .false. !< does the parameter support timeseries
  end type InputParamDefinitionType

  !> @brief Input block definition type
  !!
  !! This type is used to store information for
  !! how to read a modflow block
  !!
  !<
  type InputBlockDefinitionType
    character(len=LINELENGTH) :: blockname = '' !< name of block, e.g. DIMENSIONS
    logical(LGP) :: required = .false. !< is the block required
    logical(LGP) :: aggregate = .false. !< is this structarray style input
    logical(LGP) :: block_variable = .false. !< does this block have a block variable
    logical(LGP) :: timeseries = .false. !< does this block support timeseries
  end type InputBlockDefinitionType

end module InputDefinitionModule
