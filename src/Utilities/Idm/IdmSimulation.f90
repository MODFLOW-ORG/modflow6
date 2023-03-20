!> @brief This module contains the IdmSimulationModule
!!
!! This module contains the high-level routines for loading
!! sim namefile parameters into the input context
!!
!<
module IdmSimulationModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENMEMPATH
  use SimModule, only: store_error
  use SimVariablesModule, only: iout
  use InputOutputModule, only: openfile, getunit
  use InputDefinitionModule, only: InputParamDefinitionType
  use ModflowInputModule, only: ModflowInputType, getModflowInput
  use IdmMf6FileLoaderModule, only: input_load

  implicit none
  private
  public :: simnam_load
  public :: simnam_allocate

contains

  !> @brief MODFLOW 6 mfsim.nam input load routine
  !<
  subroutine simnam_load()
    use SimVariablesModule, only: simfile
    use GenericUtilitiesModule, only: sim_message
    integer(I4B) :: inunit
    logical :: lexist
    character(len=LINELENGTH) :: line
    !
    ! -- load mfsim.nam if it exists
    inquire (file=trim(adjustl(simfile)), exist=lexist)
    !
    if (lexist) then
      !
      ! -- write name of namfile to stdout
      write (line, '(2(1x,a))') 'Using Simulation name file:', &
        trim(adjustl(simfile))
      call sim_message(line, skipafter=1)
      !
      ! -- open namfile and load to input context
      inunit = getunit()
      call openfile(inunit, iout, trim(adjustl(simfile)), 'NAM')
      call input_load('NAM6', 'SIM', 'NAM', 'SIM', 'NAM', inunit, iout)
      !
      close (inunit)
      !
      ! -- allocate any unallocated simnam params
      call simnam_allocate()
    else
      !
      ! -- allocate  simnam params
      call simnam_allocate()
    end if
    !
    ! --return
    return
  end subroutine simnam_load

  !> @brief MODFLOW 6 mfsim.nam parameter set default value
  !<
  subroutine set_default_value(intvar, mf6varname)
    use SimVariablesModule, only: isimcontinue, isimcheck
    integer(I4B), pointer, intent(in) :: intvar
    character(len=*), intent(in) :: mf6varname
    character(len=LINELENGTH) :: errmsg
    logical(LGP) :: terminate = .true.
    !
    ! -- load defaults for keyword/integer types
    select case (mf6varname)
      !
    case ('CONTINUE')
      intvar = isimcontinue
      !
    case ('NOCHECK')
      intvar = isimcheck
      !
    case ('MAXERRORS')
      intvar = 1000 !< MessageType max_message
      !
    case ('MXITER')
      intvar = 1
      !
    case default
      write (errmsg, '(4x,a,a)') &
        '**ERROR. IdmSimulation set_default_value unhandled variable: ', &
        trim(mf6varname)
      call store_error(errmsg, terminate)
    end select
    !
    ! -- return
    return
  end subroutine set_default_value

  !> @brief MODFLOW 6 mfsim.nam input context parameter allocation
  !<
  subroutine simnam_allocate()
    use MemoryHelperModule, only: create_mem_path
    use MemoryTypeModule, only: MemoryType
    use MemoryManagerModule, only: get_isize, mem_allocate
    use SimVariablesModule, only: idm_context
    use CharacterStringModule, only: CharacterStringType
    character(len=LENMEMPATH) :: input_mempath
    type(ModflowInputType) :: mf6_input
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam, isize
    logical(LGP) :: terminate = .true.
    integer(I4B), pointer :: intvar => null()
    character(len=LINELENGTH), pointer :: cstr => null()
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: acharstr1d => null() !< variable for allocation
    character(len=LINELENGTH) :: errmsg
    !
    ! -- set memory path
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    !
    ! -- create description of input
    mf6_input = getModflowInput('NAM6', 'SIM', 'NAM', 'SIM', 'NAM')
    !
    ! -- allocate sim namfile parameters if not in input context
    do iparam = 1, size(mf6_input%param_dfns)
      !
      ! -- assign param definition pointer
      idt => mf6_input%param_dfns(iparam)
      !
      ! -- check if variable is already allocated
      call get_isize(idt%mf6varname, input_mempath, isize)
      !
      ! -- if not, allocate and set default
      if (isize < 0) then
        select case (idt%datatype)
        case ('KEYWORD', 'INTEGER')
          !
          ! -- allocate and set default
          call mem_allocate(intvar, idt%mf6varname, input_mempath)
          call set_default_value(intvar, idt%mf6varname)
          !
          ! -- reset pointer
          nullify (intvar)
        case ('STRING')
          !
          ! -- did this param originate from sim namfile RECARRAY type
          if (idt%in_record) then
            !
            ! -- allocate 0 size CharacterStringType array
            call mem_allocate(acharstr1d, LINELENGTH, 0, idt%mf6varname, &
                              input_mempath)
            !
            ! -- reset pointer
            nullify (acharstr1d)
          else
            !
            ! -- allocate empty string
            call mem_allocate(cstr, LINELENGTH, idt%mf6varname, input_mempath)
            cstr = ''
            !
            ! -- reset pointer
            nullify (cstr)
          end if
        case default
          write (errmsg, '(4x,a,a)') &
            '**ERROR. IdmSimulation unhandled datatype: ', &
            trim(idt%datatype)
          call store_error(errmsg, terminate)
        end select
      end if
    end do
    !
    ! -- return
    return
  end subroutine simnam_allocate

end module IdmSimulationModule
