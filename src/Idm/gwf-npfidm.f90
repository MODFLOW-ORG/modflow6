! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfNpfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_npf_param_definitions
  public gwf_npf_aggregate_definitions
  public gwf_npf_block_definitions
  public GwfNpfParamFoundType
  public gwf_npf_multi_package
  public gwf_npf_subpackages

  type GwfNpfParamFoundType
    logical :: ipakcb = .false.
    logical :: iprflow = .false.
    logical :: cellavg = .false.
    logical :: ithickstrt = .false.
    logical :: cvoptions = .false.
    logical :: ivarcv = .false.
    logical :: idewatcv = .false.
    logical :: iperched = .false.
    logical :: rewet_record = .false.
    logical :: irewet = .false.
    logical :: wetfct = .false.
    logical :: iwetit = .false.
    logical :: ihdwet = .false.
    logical :: xt3doptions = .false.
    logical :: ixt3d = .false.
    logical :: ixt3drhs = .false.
    logical :: ihighcellsat = .false.
    logical :: isavspdis = .false.
    logical :: isavsat = .false.
    logical :: ik22overk = .false.
    logical :: ik33overk = .false.
    logical :: tvk_filerecord = .false.
    logical :: tvk6 = .false.
    logical :: filein = .false.
    logical :: tvk6_filename = .false.
    logical :: export_ascii = .false.
    logical :: export_nc = .false.
    logical :: inewton = .false.
    logical :: satomega = .false.
    logical :: icelltype = .false.
    logical :: k = .false.
    logical :: k22 = .false.
    logical :: k33 = .false.
    logical :: angle1 = .false.
    logical :: angle2 = .false.
    logical :: angle3 = .false.
    logical :: wetdry = .false.
  end type GwfNpfParamFoundType

  logical :: gwf_npf_multi_package = .false.

  character(len=16), parameter :: &
    gwf_npf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save NPF flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print NPF flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_cellavg = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'ALTERNATIVE_CELL_AVERAGING', & ! tag name
    'CELLAVG', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'conductance weighting option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ithickstrt = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'THICKSTRT', & ! tag name
    'ITHICKSTRT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate THICKSTRT option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_cvoptions = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'CVOPTIONS', & ! tag name
    'CVOPTIONS', & ! fortran variable
    'RECORD VARIABLECV DEWATERED', & ! type
    '', & ! shape
    'vertical conductance options', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ivarcv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'VARIABLECV', & ! tag name
    'IVARCV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate VARIABLECV option', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_idewatcv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'DEWATERED', & ! tag name
    'IDEWATCV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate DEWATERED option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_iperched = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'PERCHED', & ! tag name
    'IPERCHED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate PERCHED option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_rewet_record = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'REWET_RECORD', & ! tag name
    'REWET_RECORD', & ! fortran variable
    'RECORD REWET WETFCT IWETIT IHDWET', & ! type
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
    gwfnpf_irewet = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'REWET', & ! tag name
    'IREWET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate rewetting', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_wetfct = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'WETFCT', & ! tag name
    'WETFCT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'wetting factor to use for rewetting', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_iwetit = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'IWETIT', & ! tag name
    'IWETIT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'interval to use for rewetting', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ihdwet = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'IHDWET', & ! tag name
    'IHDWET', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'flag to determine wetting equation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_xt3doptions = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'XT3DOPTIONS', & ! tag name
    'XT3DOPTIONS', & ! fortran variable
    'RECORD XT3D RHS', & ! type
    '', & ! shape
    'keyword to activate XT3D', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ixt3d = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'XT3D', & ! tag name
    'IXT3D', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate XT3D', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ixt3drhs = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'RHS', & ! tag name
    'IXT3DRHS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to XT3D on right hand side', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ihighcellsat = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'HIGHEST_CELL_SATURATION', & ! tag name
    'IHIGHCELLSAT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate HIGHEST_CELL_SATURATION option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_isavspdis = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_SPECIFIC_DISCHARGE', & ! tag name
    'ISAVSPDIS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save specific discharge', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_isavsat = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_SATURATION', & ! tag name
    'ISAVSAT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save saturation', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ik22overk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'K22OVERK', & ! tag name
    'IK22OVERK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that specified K22 is a ratio', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_ik33overk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'K33OVERK', & ! tag name
    'IK33OVERK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that specified K33 is a ratio', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_tvk_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'TVK_FILERECORD', & ! tag name
    'TVK_FILERECORD', & ! fortran variable
    'RECORD TVK6 FILEIN TVK6_FILENAME', & ! type
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
    gwfnpf_tvk6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'TVK6', & ! tag name
    'TVK6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'tvk keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
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
    gwfnpf_tvk6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'TVK6_FILENAME', & ! tag name
    'TVK6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of TVK information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_export_ascii = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_inewton = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_NEWTON', & ! tag name
    'INEWTON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'turn off Newton for unconfined cells', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_satomega = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_OMEGA', & ! tag name
    'SATOMEGA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'set saturation omega value', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_icelltype = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'ICELLTYPE', & ! tag name
    'ICELLTYPE', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'confined or convertible indicator', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_k = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'K', & ! tag name
    'K', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'hydraulic conductivity (L/T)', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_k22 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'K22', & ! tag name
    'K22', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'hydraulic conductivity of second ellipsoid axis', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_k33 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'K33', & ! tag name
    'K33', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'hydraulic conductivity of third ellipsoid axis (L/T)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_angle1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'ANGLE1', & ! tag name
    'ANGLE1', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'first anisotropy rotation angle (degrees)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_angle2 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'ANGLE2', & ! tag name
    'ANGLE2', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'second anisotropy rotation angle (degrees)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_angle3 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'ANGLE3', & ! tag name
    'ANGLE3', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'third anisotropy rotation angle (degrees)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnpf_wetdry = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NPF', & ! subcomponent
    'GRIDDATA', & ! block
    'WETDRY', & ! tag name
    'WETDRY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'wetdry threshold and factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_npf_param_definitions(*) = &
    [ &
    gwfnpf_ipakcb, &
    gwfnpf_iprflow, &
    gwfnpf_cellavg, &
    gwfnpf_ithickstrt, &
    gwfnpf_cvoptions, &
    gwfnpf_ivarcv, &
    gwfnpf_idewatcv, &
    gwfnpf_iperched, &
    gwfnpf_rewet_record, &
    gwfnpf_irewet, &
    gwfnpf_wetfct, &
    gwfnpf_iwetit, &
    gwfnpf_ihdwet, &
    gwfnpf_xt3doptions, &
    gwfnpf_ixt3d, &
    gwfnpf_ixt3drhs, &
    gwfnpf_ihighcellsat, &
    gwfnpf_isavspdis, &
    gwfnpf_isavsat, &
    gwfnpf_ik22overk, &
    gwfnpf_ik33overk, &
    gwfnpf_tvk_filerecord, &
    gwfnpf_tvk6, &
    gwfnpf_filein, &
    gwfnpf_tvk6_filename, &
    gwfnpf_export_ascii, &
    gwfnpf_export_nc, &
    gwfnpf_inewton, &
    gwfnpf_satomega, &
    gwfnpf_icelltype, &
    gwfnpf_k, &
    gwfnpf_k22, &
    gwfnpf_k33, &
    gwfnpf_angle1, &
    gwfnpf_angle2, &
    gwfnpf_angle3, &
    gwfnpf_wetdry &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_npf_aggregate_definitions(*) = &
    [ &
    InputParamDefinitionType &
    ( &
    '', & ! component
    '', & ! subcomponent
    '', & ! block
    '', & ! tag name
    '', & ! fortran variable
    '', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_npf_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwfNpfInputModule
