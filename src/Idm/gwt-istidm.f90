! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtIstInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_ist_param_definitions
  public gwt_ist_aggregate_definitions
  public gwt_ist_block_definitions
  public GwtIstParamFoundType
  public gwt_ist_multi_package
  public gwt_ist_subpackages

  type GwtIstParamFoundType
    logical :: save_flows = .false.
    logical :: budget_rec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budgetcsv_rec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: sorption = .false.
    logical :: ord1_decay = .false.
    logical :: zero_order_decay = .false.
    logical :: cim_filerecord = .false.
    logical :: cim6 = .false.
    logical :: cim6f = .false.
    logical :: cimfile = .false.
    logical :: cimprintrecord = .false.
    logical :: print_format = .false.
    logical :: formatrecord = .false.
    logical :: columns = .false.
    logical :: width = .false.
    logical :: digits = .false.
    logical :: format = .false.
    logical :: sorbate_rec = .false.
    logical :: sorbate = .false.
    logical :: sorbatefile = .false.
    logical :: export_ascii = .false.
    logical :: export_nc = .false.
    logical :: porosity = .false.
    logical :: volfrac = .false.
    logical :: zetaim = .false.
    logical :: cim = .false.
    logical :: decay = .false.
    logical :: decay_sorbed = .false.
    logical :: bulk_density = .false.
    logical :: distcoef = .false.
    logical :: sp2 = .false.
  end type GwtIstParamFoundType

  logical :: gwt_ist_multi_package = .true.

  character(len=16), parameter :: &
    gwt_ist_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtist_save_flows = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flows to budget file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_budget_rec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET_FILERECORD', & ! tag name
    'BUDGET_REC', & ! fortran variable
    'RECORD BUDGET FILEOUT BUDGETFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_budget = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET', & ! tag name
    'BUDGET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_budgetfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETFILE', & ! tag name
    'BUDGETFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_budgetcsv_rec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV_FILERECORD', & ! tag name
    'BUDGETCSV_REC', & ! fortran variable
    'RECORD BUDGETCSV FILEOUT BUDGETCSVFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_budgetcsv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV', & ! tag name
    'BUDGETCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSVFILE', & ! tag name
    'BUDGETCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_sorption = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'SORPTION', & ! tag name
    'SORPTION', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'activate sorption', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_ord1_decay = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'FIRST_ORDER_DECAY', & ! tag name
    'ORD1_DECAY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate first-order decay', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_zero_order_decay = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'ZERO_ORDER_DECAY', & ! tag name
    'ZERO_ORDER_DECAY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate zero-order decay', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_cim_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'CIM_FILERECORD', & ! tag name
    'CIM_FILERECORD', & ! fortran variable
    'RECORD CIM6 FILEOUT CIMFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_cim6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'CIM6', & ! tag name
    'CIM6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'cim6 keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_cim6f = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'CIM6F', & ! tag name
    'CIM6F', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'cim6 format keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_cimfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'CIMFILE', & ! tag name
    'CIMFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_cimprintrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'CIMPRINTRECORD', & ! tag name
    'CIMPRINTRECORD', & ! fortran variable
    'RECORD CIM6F PRINT_FORMAT FORMATRECORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_print_format = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FORMAT', & ! tag name
    'PRINT_FORMAT', & ! fortran variable
    'KEYWORD', & ! type
    ':', & ! shape
    'keyword to indicate that a print format follows', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_formatrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'FORMATRECORD', & ! tag name
    'FORMATRECORD', & ! fortran variable
    'RECORD COLUMNS WIDTH DIGITS FORMAT', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_columns = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'COLUMNS', & ! tag name
    'COLUMNS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_width = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'WIDTH', & ! tag name
    'WIDTH', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'width for each number', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_digits = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'DIGITS', & ! tag name
    'DIGITS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of digits', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_format = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'FORMAT', & ! tag name
    'FORMAT', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'write format', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_sorbate_rec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'SORBATE_FILERECORD', & ! tag name
    'SORBATE_REC', & ! fortran variable
    'RECORD SORBATE FILEOUT SORBATEFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_sorbate = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'SORBATE', & ! tag name
    'SORBATE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'sorbate keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_sorbatefile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'SORBATEFILE', & ! tag name
    'SORBATEFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_export_ascii = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_export_nc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_porosity = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'POROSITY', & ! tag name
    'POROSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'porosity of the immobile domain', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_volfrac = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'VOLFRAC', & ! tag name
    'VOLFRAC', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'volume fraction of this immobile domain', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_zetaim = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'ZETAIM', & ! tag name
    'ZETAIM', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'mass transfer rate coefficient between the mobile and immobile domains', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_cim = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'CIM', & ! tag name
    'CIM', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'initial concentration of the immobile domain', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_decay = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'DECAY', & ! tag name
    'DECAY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'first rate coefficient', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_decay_sorbed = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'DECAY_SORBED', & ! tag name
    'DECAY_SORBED', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'second rate coefficient', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_bulk_density = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'BULK_DENSITY', & ! tag name
    'BULK_DENSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'bulk density', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_distcoef = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'DISTCOEF', & ! tag name
    'DISTCOEF', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'distribution coefficient', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtist_sp2 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'IST', & ! subcomponent
    'GRIDDATA', & ! block
    'SP2', & ! tag name
    'SP2', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'second sorption parameter', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_ist_param_definitions(*) = &
    [ &
    gwtist_save_flows, &
    gwtist_budget_rec, &
    gwtist_budget, &
    gwtist_fileout, &
    gwtist_budgetfile, &
    gwtist_budgetcsv_rec, &
    gwtist_budgetcsv, &
    gwtist_budgetcsvfile, &
    gwtist_sorption, &
    gwtist_ord1_decay, &
    gwtist_zero_order_decay, &
    gwtist_cim_filerecord, &
    gwtist_cim6, &
    gwtist_cim6f, &
    gwtist_cimfile, &
    gwtist_cimprintrecord, &
    gwtist_print_format, &
    gwtist_formatrecord, &
    gwtist_columns, &
    gwtist_width, &
    gwtist_digits, &
    gwtist_format, &
    gwtist_sorbate_rec, &
    gwtist_sorbate, &
    gwtist_sorbatefile, &
    gwtist_export_ascii, &
    gwtist_export_nc, &
    gwtist_porosity, &
    gwtist_volfrac, &
    gwtist_zetaim, &
    gwtist_cim, &
    gwtist_decay, &
    gwtist_decay_sorbed, &
    gwtist_bulk_density, &
    gwtist_distcoef, &
    gwtist_sp2 &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwt_ist_aggregate_definitions(*) = &
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
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_ist_block_definitions(*) = &
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

end module GwtIstInputModule
