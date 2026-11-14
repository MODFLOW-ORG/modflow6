! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfVsetInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_vset_param_definitions
  public gwf_vset_aggregate_definitions
  public gwf_vset_block_definitions
  public GwfVsetParamFoundType
  public gwf_vset_multi_package
  public gwf_vset_is_advanced
  public gwf_vset_subpackages

  type GwfVsetParamFoundType
    logical :: et_model = .false.
    logical :: h1 = .false.
    logical :: h2 = .false.
    logical :: h3 = .false.
    logical :: h4 = .false.
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: maxbound = .false.
    logical :: cellid = .false.
    logical :: rate = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type GwfVsetParamFoundType

  logical :: gwf_vset_multi_package = .true.
  logical :: gwf_vset_is_advanced = .false.

  character(len=16), parameter :: &
    gwf_vset_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvset_et_model = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'ET_MODEL', & ! tag name
    'ET_MODEL', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model used for the root-water uptake reduction', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_h1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'H1', & ! tag name
    'H1', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'feddes anoxia pressure head', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_h2 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'H2', & ! tag name
    'H2', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'feddes upper optimal pressure head', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_h3 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'H3', & ! tag name
    'H3', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'feddes lower optimal pressure head', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_h4 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'H4', & ! tag name
    'H4', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'feddes wilting point pressure head', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'AUXMULTNAME', & ! tag name
    'AUXMULTNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of auxiliary variable for multiplier', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print transpiration rates to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save transpiration rates to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'head keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of time series information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of transpiration cells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'PERIOD', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_rate = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'PERIOD', & ! block
    'RATE', & ! tag name
    'RATE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'potential transpiration rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvset_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'transpiration name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_vset_param_definitions(*) = &
    [ &
    gwfvset_et_model, &
    gwfvset_h1, &
    gwfvset_h2, &
    gwfvset_h3, &
    gwfvset_h4, &
    gwfvset_auxiliary, &
    gwfvset_auxmultname, &
    gwfvset_boundnames, &
    gwfvset_iprpak, &
    gwfvset_iprflow, &
    gwfvset_ipakcb, &
    gwfvset_ts_filerecord, &
    gwfvset_ts6, &
    gwfvset_filein, &
    gwfvset_ts6_filename, &
    gwfvset_obs_filerecord, &
    gwfvset_obs6, &
    gwfvset_obs6_filename, &
    gwfvset_maxbound, &
    gwfvset_cellid, &
    gwfvset_rate, &
    gwfvset_auxvar, &
    gwfvset_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvset_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSET', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID RATE AUX BOUNDNAME', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_vset_aggregate_definitions(*) = &
    [ &
    gwfvset_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_vset_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfVsetInputModule
