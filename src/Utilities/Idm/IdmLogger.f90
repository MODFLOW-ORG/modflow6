!> @brief This module contains the Input Data Model Logger Module
!!
!! This module contains the subroutines for logging messages
!! to the list file as the input data model loads model input.
!!
!<
module IdmLoggerModule

  use KindModule, only: DP, LGP, I4B
  use SimVariablesModule, only: iparamlog, errmsg, idm_context
  use SimModule, only: store_error
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENCOMPONENTNAME, &
                             DISUNDEF, DIS, DISV, DISU, DIS2D, DISV1D

  implicit none
  private
  public :: idm_log_header
  public :: idm_log_close
  public :: idm_log_period_header
  public :: idm_log_period_close
  public :: idm_export
  public :: idm_log_var

  interface idm_log_var
    module procedure idm_log_var_logical, idm_log_var_int, &
      idm_log_var_int1d, idm_log_var_int2d, &
      idm_log_var_int3d, idm_log_var_dbl, &
      idm_log_var_dbl1d, idm_log_var_dbl2d, &
      idm_log_var_dbl3d, idm_log_var_str, &
      idm_log_var_ts
  end interface idm_log_var

  interface idm_export
    module procedure idm_export_int1d, idm_export_int2d, &
      idm_export_int3d, idm_export_dbl1d, &
      idm_export_dbl2d, idm_export_dbl3d
  end interface idm_export

contains

  !> @ brief log a header message
  !<
  subroutine idm_log_header(component, subcomponent, iout)
    character(len=*), intent(in) :: component !< component name
    character(len=*), intent(in) :: subcomponent !< subcomponent name
    integer(I4B), intent(in) :: iout
    if (iparamlog > 0 .and. iout > 0) then
      write (iout, '(1x,a)') 'Loading input for '//trim(component)//&
        &'/'//trim(subcomponent)
    end if
  end subroutine idm_log_header

  !> @ brief log the closing message
  !<
  subroutine idm_log_close(component, subcomponent, iout)
    character(len=*), intent(in) :: component !< component name
    character(len=*), intent(in) :: subcomponent !< subcomponent name
    integer(I4B), intent(in) :: iout
    if (iparamlog > 0 .and. iout > 0) then
      write (iout, '(1x,a)') 'Loading input complete...'
    end if
  end subroutine idm_log_close

  !> @ brief log a dynamic header message
  !<
  subroutine idm_log_period_header(component, iout)
    use TdisModule, only: kper, kstp
    character(len=*), intent(in) :: component !< component name
    integer(I4B), intent(in) :: iout
    if (iparamlog > 0 .and. iout > 0 .and. kstp == 1) then
      write (iout, '(/1x,a,i0,a)') 'IDP PERIOD ', kper, &
        ' load for component: '//trim(component)
    end if
  end subroutine idm_log_period_header

  !> @ brief log the period closing message
  !<
  subroutine idm_log_period_close(iout)
    use TdisModule, only: kstp
    integer(I4B), intent(in) :: iout
    if (iparamlog > 0 .and. iout > 0 .and. kstp == 1) then
      !backspace iout
      write (iout, '(1x,a,/)') 'IDP component dynamic load complete...'
    end if
  end subroutine idm_log_period_close

  !> @ brief log the period closing message
  !<
  subroutine idm_log_var_ts(varname, mempath, iout, is_tas)
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    logical(LGP), intent(in) :: is_tas
    if (iparamlog > 0 .and. iout > 0) then
      if (is_tas) then
        write (iout, '(3x, a, ": ", a)') &
          'Time-array-series controlled dynamic variable detected', trim(varname)
      else
        write (iout, '(3x, a, ": ", a)') &
          'Time-series controlled dynamic variable detected', trim(varname)
      end if
    end if
  end subroutine idm_log_var_ts

  !> @brief Log type specific information logical
  !<
  subroutine idm_log_var_logical(p_mem, varname, mempath, iout)
    logical(LGP), intent(in) :: p_mem !< logical scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      description = 'Logical detected'
      write (iout, '(3x, a, ": ", a, " = ", l1)') &
        trim(description), trim(varname), p_mem
    end if
  end subroutine idm_log_var_logical

  !> @brief Log type specific information integer
  !<
  subroutine idm_log_var_int(p_mem, varname, mempath, datatype, iout)
    integer(I4B), intent(in) :: p_mem !< int scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: datatype !< variable data type
    integer(I4B), intent(in) :: iout
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      if (datatype == 'KEYWORD') then
        description = 'Keyword detected'
        write (iout, '(3x, a, ": ", a)') trim(description), trim(varname)
      else
        description = 'Integer detected'
        write (iout, '(3x, a, ": ", a, " = ", i0)') &
          trim(description), trim(varname), p_mem
      end if
    end if
  end subroutine idm_log_var_int

  !> @brief Log type specific information int1d
  !<
  subroutine idm_log_var_int1d(p_mem, varname, mempath, iout)
    integer(I4B), dimension(:), contiguous, intent(in) :: p_mem !< 1d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    integer(I4B) :: min_val, max_val
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      min_val = minval(p_mem)
      max_val = maxval(p_mem)
      if (min_val == max_val) then
        description = 'Integer 1D constant array detected'
        write (iout, '(3x, a, ": ", a, " = ", i0)') &
          trim(description), trim(varname), min_val
      else
        description = 'Integer 1D array detected'
        write (iout, '(3x, a, ": ", a, a, i0, a, i0)') &
          trim(description), trim(varname), &
          ' ranges from ', min_val, ' to ', max_val
      end if
    end if
  end subroutine idm_log_var_int1d

  !> @brief Log type specific information int2d
  !<
  subroutine idm_log_var_int2d(p_mem, varname, mempath, iout)
    integer(I4B), dimension(:, :), contiguous, intent(in) :: p_mem !< 2d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    integer(I4B) :: min_val, max_val
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      min_val = minval(p_mem)
      max_val = maxval(p_mem)
      if (min_val == max_val) then
        description = 'Integer 2D constant array detected'
        write (iout, '(3x, a, ": ", a, " = ", i0)') &
          trim(description), trim(varname), min_val
      else
        description = 'Integer 2D array detected'
        write (iout, '(3x, a, ": ", a, a, i0, a, i0)') &
          trim(description), trim(varname), &
          ' ranges from ', min_val, ' to ', max_val
      end if
    end if
  end subroutine idm_log_var_int2d

  !> @brief Log type specific information int3d
  !<
  subroutine idm_log_var_int3d(p_mem, varname, mempath, iout)
    integer(I4B), dimension(:, :, :), contiguous, intent(in) :: p_mem !< 3d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    integer(I4B) :: min_val, max_val
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      min_val = minval(p_mem)
      max_val = maxval(p_mem)
      if (min_val == max_val) then
        description = 'Integer 3D constant array detected'
        write (iout, '(3x, a, ": ", a, " = ", i0)') &
          trim(description), trim(varname), min_val
      else
        description = 'Integer 3D array detected'
        write (iout, '(3x, a, ": ", a, a, i0, a, i0)') &
          trim(description), trim(varname), &
          ' ranges from ', min_val, ' to ', max_val
      end if
    end if
  end subroutine idm_log_var_int3d

  !> @brief Log type specific information double
  !<
  subroutine idm_log_var_dbl(p_mem, varname, mempath, iout)
    real(DP), intent(in) :: p_mem !< dbl scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      description = 'Double detected'
      write (iout, '(3x, a, ": ", a, " = ", G0)') &
        trim(description), trim(varname), p_mem
    end if
  end subroutine idm_log_var_dbl

  !> @brief Log type specific information dbl1d
  !<
  subroutine idm_log_var_dbl1d(p_mem, varname, mempath, iout)
    real(DP), dimension(:), contiguous, intent(in) :: p_mem !< 1d real array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    real(DP) :: min_val, max_val
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      min_val = minval(p_mem)
      max_val = maxval(p_mem)
      if (min_val == max_val) then
        description = 'Double precision 1D constant array detected'
        write (iout, '(3x, a, ": ", a, " = ", G0)') &
          trim(description), trim(varname), min_val
      else
        description = 'Double precision 1D array detected'
        write (iout, '(3x, a, ": ", a, a, G0, a, G0)') &
          trim(description), trim(varname), &
          ' ranges from ', min_val, ' to ', max_val
      end if
    end if
  end subroutine idm_log_var_dbl1d

  !> @brief Log type specific information dbl2d
  !<
  subroutine idm_log_var_dbl2d(p_mem, varname, mempath, iout)
    real(DP), dimension(:, :), contiguous, intent(in) :: p_mem !< 2d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    real(DP) :: min_val, max_val
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      min_val = minval(p_mem)
      max_val = maxval(p_mem)
      if (min_val == max_val) then
        description = 'Double precision 2D constant array detected'
        write (iout, '(3x, a, ": ", a, " = ", G0)') &
          trim(description), trim(varname), min_val
      else
        description = 'Double precision 2D array detected'
        write (iout, '(3x, a, ": ", a, a, G0, a, G0)') &
          trim(description), trim(varname), &
          ' ranges from ', min_val, ' to ', max_val
      end if
    end if
  end subroutine idm_log_var_dbl2d

  !> @brief Log type specific information dbl3d
  !<
  subroutine idm_log_var_dbl3d(p_mem, varname, mempath, iout)
    real(DP), dimension(:, :, :), contiguous, intent(in) :: p_mem !< 3d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    real(DP) :: min_val, max_val
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      min_val = minval(p_mem)
      max_val = maxval(p_mem)
      if (min_val == max_val) then
        description = 'Double precision 3D constant array detected'
        write (iout, '(3x, a, ": ", a, " = ", G0)') &
          trim(description), trim(varname), min_val
      else
        description = 'Double precision 3D array detected'
        write (iout, '(3x, a, ": ", a, a, G0, a, G0)') &
          trim(description), trim(varname), &
          ' ranges from ', min_val, ' to ', max_val
      end if
    end if
  end subroutine idm_log_var_dbl3d

  !> @brief Log type specific information str
  !<
  subroutine idm_log_var_str(p_mem, varname, mempath, iout)
    character(len=*), intent(in) :: p_mem !< pointer to str scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: iout
    character(len=LINELENGTH) :: description
    if (iparamlog > 0 .and. iout > 0) then
      description = 'String detected'
      write (iout, '(3x, a, ": ", a, " = ", a)') &
        trim(description), trim(varname), trim(p_mem)
    end if
  end subroutine idm_log_var_str

  !> @brief Create export file int1d
  !!
  !! export layered int1d parameter files
  !!
  !<
  subroutine idm_export_int1d(p_mem, varname, mempath, shapestr, iout)
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    integer(I4B), dimension(:), contiguous, intent(in) :: p_mem !< 1d integer array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: shapestr !< dfn shape string
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous :: model_shape
    integer(I4B), dimension(:, :, :), pointer, contiguous :: int3d
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    integer(I4B), dimension(3) :: dis3d_shape
    integer(I4B), dimension(2) :: dis2d_shape
    integer(I4B), pointer :: distype
    character(LENMEMPATH) :: input_mempath
    character(LENCOMPONENTNAME) :: comp, subcomp
    integer(I4B) :: i, j, k, inunit, export_dim
    logical(LGP) :: is_layered

    ! set pointer to DISENUM and MODEL_SHAPE
    call split_mem_path(mempath, comp, subcomp)
    input_mempath = create_mem_path(component=comp, context=idm_context)
    call mem_setptr(distype, 'DISENUM', input_mempath)
    call mem_setptr(model_shape, 'MODEL_SHAPE', input_mempath)

    ! set export_dim
    export_dim = distype_export_dim(distype, shapestr, is_layered)

    ! create export file(s)
    select case (export_dim)
    case (3)
      ! set reshape array
      dis3d_shape(1) = model_shape(3)
      dis3d_shape(2) = model_shape(2)
      dis3d_shape(3) = model_shape(1)
      ! allocate and reshape
      allocate (int3d(dis3d_shape(1), dis3d_shape(2), dis3d_shape(3)))
      int3d = reshape(p_mem, dis3d_shape)
      ! write export files 3D array
      do k = 1, dis3d_shape(3)
        inunit = create_export_file(varname, mempath, k, iout)
        do i = 1, model_shape(2)
          write (inunit, '(*(i0, " "))') (int3d(j, i, k), j=1, &
                                          dis3d_shape(1))
        end do
        close (inunit)
      end do
      ! cleanup
      deallocate (int3d)
    case (2)
      ! set reshape array
      dis2d_shape(1) = model_shape(2)
      dis2d_shape(2) = model_shape(1)
      ! allocate and reshape
      allocate (int2d(dis2d_shape(1), dis2d_shape(2)))
      int2d = reshape(p_mem, dis2d_shape)
      if (is_layered) then
        ! write layered export files 2D array
        do i = 1, dis2d_shape(2)
          inunit = create_export_file(varname, mempath, i, iout)
          write (inunit, '(*(i0, " "))') (int2d(j, i), j=1, dis2d_shape(1))
          close (inunit)
        end do
      else
        ! write export file 2D array
        inunit = create_export_file(varname, mempath, 0, iout)
        do i = 1, dis2d_shape(2)
          write (inunit, '(*(i0, " "))') (int2d(j, i), j=1, dis2d_shape(1))
        end do
        close (inunit)
      end if
      ! cleanup
      deallocate (int2d)
    case (1)
      ! write export file 1D array
      inunit = create_export_file(varname, mempath, 0, iout)
      write (inunit, '(*(i0, " "))') p_mem
      close (inunit)
    case default
      write (errmsg, '(a,i0)') 'EXPORT unsupported int1d export_dim=', &
        export_dim
      call store_error(errmsg, .true.)
    end select
  end subroutine idm_export_int1d

  !> @brief Create export file int2d
  !<
  subroutine idm_export_int2d(p_mem, varname, mempath, shapestr, iout)
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    integer(I4B), dimension(:, :), contiguous, intent(in) :: p_mem !< 2d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: shapestr !< dfn shape string
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous :: model_shape
    integer(I4B), pointer :: distype
    character(LENMEMPATH) :: input_mempath
    character(LENCOMPONENTNAME) :: comp, subcomp
    integer(I4B) :: i, j, inunit, export_dim
    logical(LGP) :: is_layered

    ! set pointer to DISENUM
    call split_mem_path(mempath, comp, subcomp)
    input_mempath = create_mem_path(component=comp, context=idm_context)
    call mem_setptr(distype, 'DISENUM', input_mempath)
    call mem_setptr(model_shape, 'MODEL_SHAPE', input_mempath)

    ! set export_dim
    export_dim = distype_export_dim(distype, shapestr, is_layered)

    select case (export_dim)
    case (1)
      ! write export file 1D array
      inunit = create_export_file(varname, mempath, 0, iout)
      do i = 1, size(p_mem, dim=2)
        write (inunit, '(*(i0, " "))') (p_mem(j, i), j=1, size(p_mem, dim=1))
      end do
      close (inunit)
    case (2)
      if (is_layered) then
        ! write layered export files 2D array
        do i = 1, size(p_mem, dim=2)
          inunit = create_export_file(varname, mempath, i, iout)
          write (inunit, '(*(i0, " "))') (p_mem(j, i), j=1, size(p_mem, dim=1))
          close (inunit)
        end do
      else
        ! write export file 2D array
        inunit = create_export_file(varname, mempath, 0, iout)
        do i = 1, size(p_mem, dim=2)
          write (inunit, '(*(i0, " "))') (p_mem(j, i), j=1, size(p_mem, dim=1))
        end do
        close (inunit)
      end if
    case default
      write (errmsg, '(a,i0)') 'EXPORT unsupported int2d export_dim=', &
        export_dim
      call store_error(errmsg, .true.)
    end select
  end subroutine idm_export_int2d

  !> @brief Create export file int3d
  !<
  subroutine idm_export_int3d(p_mem, varname, mempath, shapestr, iout)
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    integer(I4B), dimension(:, :, :), contiguous, intent(in) :: p_mem !< 2d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: shapestr !< dfn shape string
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous :: model_shape
    integer(I4B), pointer :: distype
    character(LENMEMPATH) :: input_mempath
    character(LENCOMPONENTNAME) :: comp, subcomp
    integer(I4B) :: i, j, k, inunit, export_dim
    logical(LGP) :: is_layered

    ! set pointer to DISENUM
    call split_mem_path(mempath, comp, subcomp)
    input_mempath = create_mem_path(component=comp, context=idm_context)
    call mem_setptr(distype, 'DISENUM', input_mempath)
    call mem_setptr(model_shape, 'MODEL_SHAPE', input_mempath)

    ! set export_dim
    export_dim = distype_export_dim(distype, shapestr, is_layered)

    select case (export_dim)
    case (3)
      ! write export files 3D array
      do k = 1, size(p_mem, dim=3)
        inunit = create_export_file(varname, mempath, k, iout)
        do i = 1, size(p_mem, dim=2)
          write (inunit, '(*(i0, " "))') (p_mem(j, i, k), j=1, size(p_mem, dim=1))
        end do
        close (inunit)
      end do
    case default
      write (errmsg, '(a,i0)') 'EXPORT unsupported int3d export_dim=', &
        export_dim
      call store_error(errmsg, .true.)
    end select
  end subroutine idm_export_int3d

  !> @brief Create export file dbl1d
  !!
  !! export layered dbl1d parameters with NODES shape
  !!
  !<
  subroutine idm_export_dbl1d(p_mem, varname, mempath, shapestr, iout)
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    real(DP), dimension(:), contiguous, intent(in) :: p_mem !< 1d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: shapestr !< dfn shape string
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous :: model_shape
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B), dimension(3) :: dis3d_shape
    integer(I4B), dimension(2) :: dis2d_shape
    integer(I4B), pointer :: distype
    character(LENMEMPATH) :: input_mempath
    character(LENCOMPONENTNAME) :: comp, subcomp
    integer(I4B) :: i, j, k, inunit, export_dim
    logical(LGP) :: is_layered

    ! set pointer to DISENUM and MODEL_SHAPE
    call split_mem_path(mempath, comp, subcomp)
    input_mempath = create_mem_path(component=comp, context=idm_context)
    call mem_setptr(distype, 'DISENUM', input_mempath)
    call mem_setptr(model_shape, 'MODEL_SHAPE', input_mempath)

    ! set export_dim
    export_dim = distype_export_dim(distype, shapestr, is_layered)

    ! create export file(s)
    select case (export_dim)
    case (3)
      ! set reshape array
      dis3d_shape(1) = model_shape(3)
      dis3d_shape(2) = model_shape(2)
      ! allocate and reshape
      dis3d_shape(3) = model_shape(1)
      allocate (dbl3d(dis3d_shape(1), dis3d_shape(2), dis3d_shape(3)))
      dbl3d = reshape(p_mem, dis3d_shape)
      do k = 1, dis3d_shape(3)
        ! write export files 3D array
        inunit = create_export_file(varname, mempath, k, iout)
        do i = 1, model_shape(2)
          write (inunit, '(*(G0.10, " "))') (dbl3d(j, i, k), j=1, &
                                             dis3d_shape(1))
        end do
        close (inunit)
      end do
      ! cleanup
      deallocate (dbl3d)
    case (2)
      ! set reshape array
      dis2d_shape(1) = model_shape(2)
      dis2d_shape(2) = model_shape(1)
      ! allocate and reshape
      allocate (dbl2d(dis2d_shape(1), dis2d_shape(2)))
      dbl2d = reshape(p_mem, dis2d_shape)
      if (is_layered) then
        ! write layered export files 2D array
        do i = 1, dis2d_shape(2)
          inunit = create_export_file(varname, mempath, i, iout)
          write (inunit, '(*(G0.10, " "))') (dbl2d(j, i), j=1, dis2d_shape(1))
          close (inunit)
        end do
      else
        ! write export file 2D array
        inunit = create_export_file(varname, mempath, 0, iout)
        do i = 1, dis2d_shape(2)
          write (inunit, '(*(G0.10, " "))') (dbl2d(j, i), j=1, dis2d_shape(1))
        end do
        close (inunit)
      end if
      ! cleanup
      deallocate (dbl2d)
    case (1)
      ! write export file 1D array
      inunit = create_export_file(varname, mempath, 0, iout)
      write (inunit, '(*(G0.10, " "))') p_mem
      close (inunit)
    case default
      write (errmsg, '(a,i0)') 'EXPORT unsupported dbl1d export_dim=', &
        export_dim
      call store_error(errmsg, .true.)
    end select
  end subroutine idm_export_dbl1d

  !> @brief Create export file dbl2d
  !<
  subroutine idm_export_dbl2d(p_mem, varname, mempath, shapestr, iout)
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    real(DP), dimension(:, :), contiguous, intent(in) :: p_mem !< 2d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: shapestr !< dfn shape string
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous :: model_shape
    integer(I4B), pointer :: distype
    character(LENMEMPATH) :: input_mempath
    character(LENCOMPONENTNAME) :: comp, subcomp
    integer(I4B) :: i, j, inunit, export_dim
    logical(LGP) :: is_layered

    ! set pointer to DISENUM
    call split_mem_path(mempath, comp, subcomp)
    input_mempath = create_mem_path(component=comp, context=idm_context)
    call mem_setptr(distype, 'DISENUM', input_mempath)
    call mem_setptr(model_shape, 'MODEL_SHAPE', input_mempath)

    ! set export_dim
    export_dim = distype_export_dim(distype, shapestr, is_layered)

    select case (export_dim)
    case (1)
      ! write export file 1D array
      inunit = create_export_file(varname, mempath, 0, iout)
      do i = 1, size(p_mem, dim=2)
        write (inunit, '(*(G0.10, " "))') (p_mem(j, i), j=1, size(p_mem, dim=1))
      end do
      close (inunit)
    case (2)
      if (is_layered) then
        ! write layered export files 2D array
        do i = 1, size(p_mem, dim=2)
          inunit = create_export_file(varname, mempath, i, iout)
          write (inunit, '(*(G0.10, " "))') (p_mem(j, i), j=1, size(p_mem, dim=1))
          close (inunit)
        end do
      else
        ! write export file 2D array
        inunit = create_export_file(varname, mempath, 0, iout)
        do i = 1, size(p_mem, dim=2)
          write (inunit, '(*(G0.10, " "))') (p_mem(j, i), j=1, size(p_mem, dim=1))
        end do
        close (inunit)
      end if
    case default
      write (errmsg, '(a,i0)') 'EXPORT unsupported dbl2d export_dim=', &
        export_dim
      call store_error(errmsg, .true.)
    end select
  end subroutine idm_export_dbl2d

  !> @brief Create export file dbl3d
  !<
  subroutine idm_export_dbl3d(p_mem, varname, mempath, shapestr, iout)
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    real(DP), dimension(:, :, :), contiguous, intent(in) :: p_mem !< 2d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    character(len=*), intent(in) :: shapestr !< dfn shape string
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous :: model_shape
    integer(I4B), pointer :: distype
    character(LENMEMPATH) :: input_mempath
    character(LENCOMPONENTNAME) :: comp, subcomp
    integer(I4B) :: i, j, k, inunit, export_dim
    logical(LGP) :: is_layered

    ! set pointer to DISENUM
    call split_mem_path(mempath, comp, subcomp)
    input_mempath = create_mem_path(component=comp, context=idm_context)
    call mem_setptr(distype, 'DISENUM', input_mempath)
    call mem_setptr(model_shape, 'MODEL_SHAPE', input_mempath)

    ! set export_dim
    export_dim = distype_export_dim(distype, shapestr, is_layered)

    select case (export_dim)
    case (3)
      ! write export files 3D array
      do k = 1, size(p_mem, dim=3)
        inunit = create_export_file(varname, mempath, k, iout)
        do i = 1, size(p_mem, dim=2)
          write (inunit, '(*(G0.10, " "))') (p_mem(j, i, k), j=1, &
                                             size(p_mem, dim=1))
        end do
        close (inunit)
      end do
    case default
      write (errmsg, '(a,i0)') 'EXPORT unsupported dbl3d export_dim=', &
        export_dim
      call store_error(errmsg, .true.)
    end select
  end subroutine idm_export_dbl3d

  !> @brief Set dis type export_dim
  !!
  !! Set the dimension of the export
  !<
  function distype_export_dim(distype, shapestr, is_layered) &
    result(export_dim)
    integer(I4B), pointer, intent(in) :: distype
    character(len=*), intent(in) :: shapestr !< dfn shape string
    logical(LGP), intent(inout) :: is_layered !< does this data represent layers
    integer(I4B) :: export_dim

    ! initialize is_layered to false
    is_layered = .false.

    select case (distype)
    case (DIS)
      if (shapestr == 'NODES') then
        export_dim = 3
        is_layered = .true.
      else if (shapestr == 'NCOL NROW NLAY') then
        export_dim = 3
        is_layered = .true.
      else
        export_dim = 1
      end if
    case (DISV)
      if (shapestr == 'NODES') then
        export_dim = 2
        is_layered = .true.
      else if (shapestr == 'NCPL NLAY') then
        export_dim = 2
        is_layered = .true.
      else
        export_dim = 1
      end if
    case (DIS2D)
      if (shapestr == 'NODES') then
        export_dim = 2
      else if (shapestr == 'NCOL NROW') then
        export_dim = 2
      else
        export_dim = 1
      end if
    case (DISU, DISV1D)
      export_dim = 1
    case default
      export_dim = 0
    end select
  end function distype_export_dim

  !> @brief Create export file
  !!
  !! Name formats where l=layer, a=auxiliary, p=period
  !!    : <comp>-<subcomp>.varname.txt
  !!    : <comp>-<subcomp>.varname.l<num>.txt
  !!    : <comp>-<subcomp>.varname.p<num>.txt
  !!    : <comp>-<subcomp>.varname.a<num>.p<num>.txt
  !<
  function create_export_file(varname, mempath, layer, iout) &
    result(inunit)
    use ConstantsModule, only: LENVARNAME
    use InputOutputModule, only: openfile, getunit
    use InputOutputModule, only: upcase, lowcase
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable memory path
    integer(I4B), intent(in) :: layer
    integer(I4B), intent(in) :: iout
    integer(I4B) :: inunit
    character(len=LENCOMPONENTNAME) :: comp, subcomp
    character(len=LINELENGTH) :: filename, suffix

    ! split the mempath
    call split_mem_path(mempath, comp, subcomp)
    call lowcase(comp)
    call lowcase(subcomp)

    ! build suffix
    suffix = varname
    call lowcase(suffix)
    if (layer > 0) then
      write (suffix, '(a,i0)') trim(suffix)//'.l', layer
    end if
    suffix = trim(suffix)//'.txt'

    ! set filename
    filename = trim(comp)//'-'//trim(subcomp)//'.'//trim(suffix)

    ! silently create the array file
    inunit = getunit()
    call openfile(inunit, 0, filename, 'EXPORT', filstat_opt='REPLACE')
  end function create_export_file

end module IdmLoggerModule
