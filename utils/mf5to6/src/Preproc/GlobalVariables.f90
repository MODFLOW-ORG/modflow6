module GlobalVariablesModule
  
  use ConstantsModule, only: MAXCHARLEN
  use VersionModule, only: FULLVERSION
  use OpenSpecModule,  only: ACCESS, ACTION, FORM

  implicit none

  private
  public :: prognamconv, prognamlong, mfvnam, ilgr, ilunit, ngrids, &
            NIUNIT, cunit, verbose, LgrBilinear, masteridomain, &
            GetNextIgrid, optfile, PathToPostObsMf, ScriptType, echo, &
            msgc

  character(len=60) :: prognamconv, prognamlong
  character(len=48) :: mfvnam
  parameter (prognamconv='Mf5to6')
  parameter (prognamlong=trim(prognamconv)//' - Converter for MODFLOW (2005, NWT, LGR) to MODFLOW 6')
  parameter (mfvnam='Version ' // FULLVERSION)
  integer  :: ngrid = 0
  integer  :: ilgr, ilunit, ngrids
  integer, parameter :: NIUNIT=100
  character(len=4) :: cunit(NIUNIT)
  character(len=MAXCHARLEN) :: optfile='', PathToPostObsMf=''
  character(len=MAXCHARLEN) :: msgc = 'Conversion successful!'
  character(len=10) :: ScriptType='BATCH'
  logical :: echo = .false.
  logical :: verbose = .false.
  logical :: LgrBilinear = .false.
  integer, dimension(:,:,:), pointer :: masteridomain => null()
  data cunit/'BCF6', 'WEL ', 'DRN ', 'RIV ', 'EVT ', 'gfd ', 'GHB ', & !  7
             'RCH ', 'SIP ', 'DE4 ', '    ', 'OC  ', 'PCG ', 'lmg ', & ! 14
             'gwt ', 'FHB ', 'RES ', 'STR ', 'IBS ', 'CHD ', 'HFB6', & ! 21
             'LAK ', 'LPF ', 'DIS ', '    ', 'PVAL', '    ', 'HOB ', & ! 28
             '    ', '    ', 'ZONE', 'MULT', 'DROB', 'RVOB', 'GBOB', & ! 35
             'STOB', 'HUF2', 'CHOB', 'ETS ', 'DRT ', '    ', 'GMG ', & ! 42
             'HYD ', 'SFR ', '    ', 'GAGE', 'LVDA', '    ', 'LMT6', & ! 49
             'MNW2', 'MNWI', 'MNW1', 'KDEP', 'SUB ', 'UZF ', 'gwm ', & ! 56
             'SWT ', 'cfp ', 'PCGN', '    ', '    ', 'UPW ', 'NWT ', & ! 63
             '    ', 'SWI2', '    ', '    ', '    ', '    ', '    ', & ! 70
             30*'    '/

contains

  integer function GetNextIgrid()
    implicit none
    !
    ngrid = ngrid + 1
    GetNextIgrid = ngrid
    !
    return
  end function GetNextIgrid

end module GlobalVariablesModule
