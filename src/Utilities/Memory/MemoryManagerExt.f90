module MemoryManagerExtModule

  use KindModule, only: DP, LGP, I4B, I8B
  use SimModule, only: store_error
  use ConstantsModule, only: MVALIDATE
  use SimVariablesModule, only: isim_mode
  use MemoryTypeModule, only: MemoryType
  use MemoryManagerModule, only: memorystore, get_from_memorystore, mem_release
  use MemoryContainerIteratorModule, only: MemoryContainerIteratorType

  implicit none
  private
  public :: mem_set_value
  public :: memorystore_remove
  public :: memorystore_release

  interface mem_set_value
    module procedure mem_set_value_logical, mem_set_value_int, &
      mem_set_value_int_setval, mem_set_value_str_mapped_int, &
      mem_set_value_logical1d, mem_set_value_logical1d_mapped, &
      mem_set_value_int1d, mem_set_value_int1d_mapped, &
      mem_set_value_int2d, mem_set_value_int3d, mem_set_value_dbl, &
      mem_set_value_dbl1d, mem_set_value_dbl1d_mapped, &
      mem_set_value_dbl2d, mem_set_value_dbl3d, mem_set_value_str, &
      mem_set_value_charstr1d
  end interface mem_set_value

contains

  subroutine memorystore_remove(component, subcomponent, context)
    use MemoryHelperModule, only: create_mem_path
    use ConstantsModule, only: LENMEMPATH
    character(len=*), intent(in) :: component !< name of the solution, model, or exchange
    character(len=*), intent(in), optional :: subcomponent !< name of the package (optional)
    character(len=*), intent(in), optional :: context !< name of the context (optional)
    character(len=LENMEMPATH) :: memory_path !< the memory path
    type(MemoryType), pointer :: mt
    type(MemoryContainerIteratorType), allocatable :: itr
    logical(LGP) :: removed

    memory_path = create_mem_path(component, subcomponent, context)
    removed = .true. !< initialize the loop

    do while (removed)
      removed = .false.
      itr = memorystore%iterator()
      do while (itr%has_next())
        call itr%next()
        mt => itr%value()
        ! guard: mt_associated() is false for entries already released via
        ! mem_set_value(..., release=.true.) prior to this call
        if (mt%path == memory_path .and. mt%mt_associated()) then
          call mt%mt_deallocate()
          removed = .true.
          deallocate (itr)
          exit
        end if
      end do
    end do
  end subroutine memorystore_remove

  !> @brief Release a single variable from the memory store
  !!
  !! Looks up the variable by name and path and deallocates its data.
  !! Safe to call when the variable is not found or was already released.
  !<
  subroutine memorystore_release(varname, memory_path)
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    logical(LGP) :: checkfail = .false.

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: already released
    if (isim_mode /= MVALIDATE) call mem_release(mt)
  end subroutine memorystore_release

  !> @brief Set pointer to value of memory list logical variable
  !<
  subroutine mem_set_value_logical(p_mem, varname, memory_path, found, release)
    logical(LGP), pointer, intent(inout) :: p_mem !< pointer to logical scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'INTEGER') then
      if (mt%intsclr == 0) then
        p_mem = .false.
      else
        p_mem = .true.
      end if
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_logical

  !> @brief Set pointer to value of memory list int variable
  !<
  subroutine mem_set_value_int(p_mem, varname, memory_path, found, release)
    integer(I4B), pointer, intent(inout) :: p_mem !< pointer to int scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'INTEGER') then
      p_mem = mt%intsclr
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_int

  subroutine mem_set_value_int_setval(p_mem, varname, memory_path, setval, &
                                      found, release)
    integer(I4B), pointer, intent(inout) :: p_mem !< pointer to int scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    integer(I4B), intent(in) :: setval !< set p_mem to setval if varname found
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released

    p_mem = setval
    if (do_release) call mem_release(mt)
  end subroutine mem_set_value_int_setval

  subroutine mem_set_value_str_mapped_int(p_mem, varname, memory_path, &
                                          str_list, found, release)
    integer(I4B), pointer, intent(inout) :: p_mem !< pointer to int scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    character(len=*), dimension(:), intent(in) :: str_list
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: i

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'STRING') then
      do i = 1, size(str_list)
        if (mt%strsclr == str_list(i)) then
          p_mem = i
        end if
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_str_mapped_int

  !> @brief Set pointer to value of memory list 1d logical array variable
  !<
  subroutine mem_set_value_logical1d(p_mem, varname, memory_path, found, &
                                     release)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: p_mem !< pointer to 1d logical array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'LOGICAL') then
      if (size(mt%alogical1d) /= size(p_mem)) then
        call store_error('mem_set_value() size mismatch logical1d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do n = 1, size(mt%alogical1d)
        p_mem(n) = mt%alogical1d(n)
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_logical1d

  !> @brief Set pointer to value of memory list 1d logical array variable with mapping
  !<
  subroutine mem_set_value_logical1d_mapped(p_mem, varname, memory_path, map, &
                                            found, release)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: p_mem !< pointer to 1d logical array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: map !< pointer to 1d int mapping array
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'LOGICAL') then
      if (associated(map)) then
        do n = 1, size(p_mem)
          p_mem(n) = mt%alogical1d(map(n))
        end do
      else
        if (size(mt%alogical1d) /= size(p_mem)) then
          call store_error('mem_set_value() size mismatch logical1d, varname='//&
                           &trim(varname), terminate=.TRUE.)
        end if
        do n = 1, size(mt%alogical1d)
          p_mem(n) = mt%alogical1d(n)
        end do
      end if
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_logical1d_mapped

  !> @brief Set pointer to value of memory list 1d int array variable
  !<
  subroutine mem_set_value_int1d(p_mem, varname, memory_path, found, release)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: p_mem !< pointer to 1d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'INTEGER') then
      if (size(mt%aint1d) /= size(p_mem)) then
        call store_error('mem_set_value() size mismatch int1d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do n = 1, size(mt%aint1d)
        p_mem(n) = mt%aint1d(n)
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_int1d

  !> @brief Set pointer to value of memory list 1d int array variable with mapping
  !<
  subroutine mem_set_value_int1d_mapped(p_mem, varname, memory_path, map, &
                                        found, release)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: p_mem !< pointer to 1d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: map !< pointer to 1d int mapping array
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'INTEGER') then
      if (associated(map)) then
        do n = 1, size(p_mem)
          p_mem(n) = mt%aint1d(map(n))
        end do
      else
        if (size(mt%aint1d) /= size(p_mem)) then
          call store_error('mem_set_value() size mismatch int1d, varname='//&
                           &trim(varname), terminate=.TRUE.)
        end if
        do n = 1, size(mt%aint1d)
          p_mem(n) = mt%aint1d(n)
        end do
      end if
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_int1d_mapped

  !> @brief Set pointer to value of memory list 2d int array variable
  !<
  subroutine mem_set_value_int2d(p_mem, varname, memory_path, found, release)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: p_mem !< pointer to 2d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: i, j

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'INTEGER') then
      if (size(mt%aint2d, dim=1) /= size(p_mem, dim=1) .or. &
          size(mt%aint2d, dim=2) /= size(p_mem, dim=2)) then
        call store_error('mem_set_value() size mismatch int2d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do j = 1, size(mt%aint2d, dim=2)
        do i = 1, size(mt%aint2d, dim=1)
          p_mem(i, j) = mt%aint2d(i, j)
        end do
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_int2d

  !> @brief Set pointer to value of memory list 3d int array variable
  !<
  subroutine mem_set_value_int3d(p_mem, varname, memory_path, found, release)
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(inout) :: p_mem !< pointer to 3d int array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: i, j, k

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'INTEGER') then
      if (size(mt%aint3d, dim=1) /= size(p_mem, dim=1) .or. &
          size(mt%aint3d, dim=2) /= size(p_mem, dim=2) .or. &
          size(mt%aint3d, dim=3) /= size(p_mem, dim=3)) then
        call store_error('mem_set_value() size mismatch int3d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do k = 1, size(mt%aint3d, dim=3)
        do j = 1, size(mt%aint3d, dim=2)
          do i = 1, size(mt%aint3d, dim=1)
            p_mem(i, j, k) = mt%aint3d(i, j, k)
          end do
        end do
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_int3d

  !> @brief Set pointer to value of memory list double variable
  !<
  subroutine mem_set_value_dbl(p_mem, varname, memory_path, found, release)
    real(DP), pointer, intent(inout) :: p_mem !< pointer to dbl scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'DOUBLE') then
      p_mem = mt%dblsclr
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_dbl

  !> @brief Set pointer to value of memory list 1d dbl array variable
  !<
  subroutine mem_set_value_dbl1d(p_mem, varname, memory_path, found, release)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: p_mem !< pointer to 1d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'DOUBLE') then
      if (size(mt%adbl1d) /= size(p_mem)) then
        call store_error('mem_set_value() size mismatch dbl1d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do n = 1, size(mt%adbl1d)
        p_mem(n) = mt%adbl1d(n)
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_dbl1d

  !> @brief Set pointer to value of memory list 1d dbl array variable with mapping
  !<
  subroutine mem_set_value_dbl1d_mapped(p_mem, varname, memory_path, map, &
                                        found, release)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: p_mem !< pointer to 1d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: map !< pointer to 1d int mapping array
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'DOUBLE') then
      if (associated(map)) then
        do n = 1, size(p_mem)
          p_mem(n) = mt%adbl1d(map(n))
        end do
      else
        if (size(mt%adbl1d) /= size(p_mem)) then
          call store_error('mem_set_value() size mismatch dbl1d, varname='//&
                           &trim(varname), terminate=.TRUE.)
        end if
        do n = 1, size(mt%adbl1d)
          p_mem(n) = mt%adbl1d(n)
        end do
      end if
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_dbl1d_mapped

  !> @brief Set pointer to value of memory list 2d dbl array variable
  !<
  subroutine mem_set_value_dbl2d(p_mem, varname, memory_path, found, release)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: p_mem !< pointer to 2d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: i, j

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'DOUBLE') then
      if (size(mt%adbl2d, dim=1) /= size(p_mem, dim=1) .or. &
          size(mt%adbl2d, dim=2) /= size(p_mem, dim=2)) then
        call store_error('mem_set_value() size mismatch dbl2d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do j = 1, size(mt%adbl2d, dim=2)
        do i = 1, size(mt%adbl2d, dim=1)
          p_mem(i, j) = mt%adbl2d(i, j)
        end do
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_dbl2d

  !> @brief Set pointer to value of memory list 3d dbl array variable
  !<
  subroutine mem_set_value_dbl3d(p_mem, varname, memory_path, found, release)
    real(DP), dimension(:, :, :), pointer, contiguous, intent(inout) :: p_mem !< pointer to 3d dbl array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: i, j, k

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'DOUBLE') then
      if (size(mt%adbl3d, dim=1) /= size(p_mem, dim=1) .or. &
          size(mt%adbl3d, dim=2) /= size(p_mem, dim=2) .or. &
          size(mt%adbl3d, dim=3) /= size(p_mem, dim=3)) then
        call store_error('mem_set_value() size mismatch dbl3d, varname='//&
                         &trim(varname), terminate=.TRUE.)
      end if
      do k = 1, size(mt%adbl3d, dim=3)
        do j = 1, size(mt%adbl3d, dim=2)
          do i = 1, size(mt%adbl3d, dim=1)
            p_mem(i, j, k) = mt%adbl3d(i, j, k)
          end do
        end do
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_dbl3d

  subroutine mem_set_value_str(p_mem, varname, memory_path, found, release)
    character(len=*), intent(inout) :: p_mem !< pointer to str scalar
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'STRING') then
      p_mem = mt%strsclr
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_str

  subroutine mem_set_value_charstr1d(p_mem, varname, memory_path, found, &
                                     release)
    use CharacterStringModule, only: CharacterStringType
    type(CharacterStringType), dimension(:), &
      pointer, contiguous, intent(inout) :: p_mem !< pointer to charstr 1d array
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: memory_path !< path where variable is stored
    logical(LGP), intent(inout) :: found
    logical(LGP), intent(in), optional :: release !< if true (default), deallocate input context memory after copy
    type(MemoryType), pointer :: mt
    logical(LGP) :: checkfail = .false.
    logical(LGP) :: do_release
    integer(I4B) :: n

    do_release = (isim_mode /= MVALIDATE)
    if (present(release)) do_release = release

    call get_from_memorystore(varname, memory_path, mt, found, checkfail)
    if (.not. found) return
    if (.not. mt%mt_associated()) return ! guard: entry was previously released
    if (mt%memtype(1:index(mt%memtype, ' ')) == 'STRING') then
      do n = 1, size(mt%acharstr1d)
        p_mem(n) = mt%acharstr1d(n)
      end do
      if (do_release) call mem_release(mt)
    end if
  end subroutine mem_set_value_charstr1d

end module MemoryManagerExtModule
