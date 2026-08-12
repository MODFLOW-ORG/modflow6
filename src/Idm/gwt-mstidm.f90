! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtMstInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_mst_param_definitions
  public gwt_mst_aggregate_definitions
  public gwt_mst_block_definitions
  public GwtMstParamFoundType
  public gwt_mst_multi_package
  public gwt_mst_subpackages

  type GwtMstParamFoundType
    logical :: save_flows = .false.
    logical :: order1_decay = .false.
    logical :: order0_decay = .false.
    logical :: sorption = .false.
    logical :: istrand = .false.
    logical :: stranded_rec = .false.
    logical :: stranded = .false.
    logical :: strandedfile = .false.
    logical :: sorbate_rec = .false.
    logical :: sorbate = .false.
    logical :: fileout = .false.
    logical :: sorbatefile = .false.
    logical :: export_ascii = .false.
    logical :: export_nc = .false.
    logical :: porosity = .false.
    logical :: decay = .false.
    logical :: decay_sorbed = .false.
    logical :: bulk_density = .false.
    logical :: distcoef = .false.
    logical :: sp2 = .false.
  end type GwtMstParamFoundType

  logical :: gwt_mst_multi_package = .false.

  character(len=16), parameter :: &
    gwt_mst_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtmst_save_flows = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_order1_decay = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'FIRST_ORDER_DECAY', & ! tag name
    'ORDER1_DECAY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate first-order decay', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_order0_decay = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'ZERO_ORDER_DECAY', & ! tag name
    'ORDER0_DECAY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate zero-order decay', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_sorption = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'SORPTION', & ! tag name
    'SORPTION', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'activate sorption', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_istrand = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'STRANDED_MASS', & ! tag name
    'ISTRAND', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate stranded mass', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_stranded_rec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'STRANDED_FILERECORD', & ! tag name
    'STRANDED_REC', & ! fortran variable
    'RECORD STRANDED FILEOUT STRANDEDFILE', & ! type
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
    gwtmst_stranded = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'STRANDED', & ! tag name
    'STRANDED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'stranded keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_strandedfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'STRANDEDFILE', & ! tag name
    'STRANDEDFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_sorbate_rec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'SORBATE_FILERECORD', & ! tag name
    'SORBATE_REC', & ! fortran variable
    'RECORD SORBATE FILEOUT SORBATEFILE', & ! type
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
    gwtmst_sorbate = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'SORBATE', & ! tag name
    'SORBATE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'sorbate keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
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
    gwtmst_sorbatefile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'OPTIONS', & ! block
    'SORBATEFILE', & ! tag name
    'SORBATEFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_export_ascii = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
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
    gwtmst_export_nc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
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
    gwtmst_porosity = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'GRIDDATA', & ! block
    'POROSITY', & ! tag name
    'POROSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'porosity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_decay = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'GRIDDATA', & ! block
    'DECAY', & ! tag name
    'DECAY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'aqueous phase decay rate coefficient', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_decay_sorbed = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'GRIDDATA', & ! block
    'DECAY_SORBED', & ! tag name
    'DECAY_SORBED', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'sorbed phase decay rate coefficient', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_bulk_density = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'GRIDDATA', & ! block
    'BULK_DENSITY', & ! tag name
    'BULK_DENSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'bulk density', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_distcoef = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'GRIDDATA', & ! block
    'DISTCOEF', & ! tag name
    'DISTCOEF', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'distribution coefficient', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmst_sp2 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MST', & ! subcomponent
    'GRIDDATA', & ! block
    'SP2', & ! tag name
    'SP2', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'second sorption parameter', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_mst_param_definitions(*) = &
    [ &
    gwtmst_save_flows, &
    gwtmst_order1_decay, &
    gwtmst_order0_decay, &
    gwtmst_sorption, &
    gwtmst_istrand, &
    gwtmst_stranded_rec, &
    gwtmst_stranded, &
    gwtmst_strandedfile, &
    gwtmst_sorbate_rec, &
    gwtmst_sorbate, &
    gwtmst_fileout, &
    gwtmst_sorbatefile, &
    gwtmst_export_ascii, &
    gwtmst_export_nc, &
    gwtmst_porosity, &
    gwtmst_decay, &
    gwtmst_decay_sorbed, &
    gwtmst_bulk_density, &
    gwtmst_distcoef, &
    gwtmst_sp2 &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwt_mst_aggregate_definitions(*) = &
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
    gwt_mst_block_definitions(*) = &
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

end module GwtMstInputModule
