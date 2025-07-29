module GwfExgMoverModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LENMODELNAME, LENPACKAGENAME, DZERO, DNODATA
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryHelperModule, only: split_mem_path
  use VirtualModelModule
  use BaseDisModule
  use GwfMvrModule
  use PackageMoverModule, only: PackageMoverType, set_packagemover_pointer
  implicit none
  private

  public :: exg_mvr_cr

  !> @brief Extends model mover for exchanges to also handle the
  !< parallel case where the models are not on the same process.
  type, public, extends(GwfMvrType) :: GwfExgMoverType
    class(VirtualModelType), pointer :: model1 => null() !< virtual model 1
    class(VirtualModelType), pointer :: model2 => null() !< virtual model 2
    logical(LGP), dimension(:), pointer, contiguous :: prov_is_m1 => null() !< .true. when the providing package is part of model 1
    real(DP), dimension(:), pointer, contiguous :: qpactual_m1 => null() !< stores qpactual for synchronization when provider is in model 1
    real(DP), dimension(:), pointer, contiguous :: qpactual_m2 => null() !< stores qpactual for synchronization when provider is in model 2
    real(DP), dimension(:), pointer, contiguous :: qavailable_m1 => null() !< stores qavailable for synchronization when provider is in model 1
    real(DP), dimension(:), pointer, contiguous :: qavailable_m2 => null() !< stores qavailable for synchronization when provider is in model 2
    integer(I4B), dimension(:), pointer, contiguous :: id_mapped_m1 => null() !< stores the mapped feature ids for synchronization when provider is in model 1
    integer(I4B), dimension(:), pointer, contiguous :: id_mapped_m2 => null() !< stores the mapped feature ids for synchronization when provider is in model 2
  contains
    procedure :: mvr_da => xmvr_da
    procedure :: xmvr_cf
    procedure :: mvr_fc => xmvr_fc
    procedure :: mvr_bd => xmvr_bd
    procedure :: check_packages => xmvr_check_packages
    procedure :: assign_packagemovers => xmvr_assign_packagemovers
    procedure :: initialize_movers => xmvr_initialize_movers
    procedure :: allocate_arrays => xmvr_allocate_arrays
  end type

contains

  subroutine exg_mvr_cr(exg_mvr, name_parent, input_mempath, inunit, iout, dis)
    class(GwfExgMoverType), pointer :: exg_mvr
    character(len=*), intent(in) :: name_parent
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer :: dis

    allocate (exg_mvr)

    ! Init through base
    call exg_mvr%mvr_init(name_parent, input_mempath, inunit, iout, dis, 1)

  end subroutine exg_mvr_cr

  subroutine xmvr_check_packages(this)
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_setptr
    use SimModule, only: store_error, count_errors, store_error_unit
    class(GwfExgMoverType), intent(inout) :: this
    ! local
    character(len=LENMODELNAME) :: mname
    character(len=LENPACKAGENAME) :: pname
    class(VirtualModelType), pointer :: vm
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: i
    integer(I4B), pointer :: imover_ptr

    do i = 1, size(this%pckMemPaths)
      ! check only when local
      call split_mem_path(this%pckMemPaths(i), mname, pname)
      vm => get_virtual_model(mname)
      if (vm%is_local) then
        ! check if PackageMover is active in package:
        imover_ptr => null()
        call mem_setptr(imover_ptr, 'IMOVER', trim(this%pckMemPaths(i)))
        if (imover_ptr == 0) then
          write (errmsg, '(a, a, a)') &
            'ERROR.  MODEL AND PACKAGE "', &
            trim(this%pckMemPaths(i)), &
            '" DOES NOT HAVE MOVER SPECIFIED IN OPTIONS BLOCK.'
          call store_error(errmsg)
        end if
      end if
    end do

    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if

  end subroutine xmvr_check_packages

  !> @brief Overrides GWF MVR routine to skip assigning
  !< pointers when the package is not local
  subroutine xmvr_assign_packagemovers(this)
    class(GwfExgMoverType), intent(inout) :: this !< this exchange mover
    ! local
    integer(I4B) :: i
    character(len=LENMODELNAME) :: mname
    character(len=LENPACKAGENAME) :: pname
    class(VirtualModelType), pointer :: vm

    do i = 1, size(this%pckMemPaths)
      if (this%pakmovers(i)%memoryPath == '') then
        ! is it local?
        call split_mem_path(this%pckMemPaths(i), mname, pname)
        vm => get_virtual_model(mname)
        if (vm%is_local) then
          ! yes, we need the pointers
          call set_packagemover_pointer(this%pakmovers(i), &
                                        trim(this%pckMemPaths(i)))
        end if
      end if
    end do
  end subroutine xmvr_assign_packagemovers

  !> @brief Overrides mover initialization in GWF MVR to
  !! deactivate remote parts and build up sync. arrays
  !< for mapped feature ids
  subroutine xmvr_initialize_movers(this, nr_active_movers)
    class(GwfExgMoverType) :: this
    integer(I4B) :: nr_active_movers
    ! local
    integer(I4B) :: i
    character(len=LENMODELNAME) :: mname
    character(len=LENPACKAGENAME) :: pname
    class(VirtualModelType), pointer :: vm
    class(PackageMoverType), allocatable :: pkg_mvr

    call this%GwfMvrType%initialize_movers(nr_active_movers)

    this%prov_is_m1 = .false.

    ! deactivate remote parts
    do i = 1, nr_active_movers
      call split_mem_path(this%mvr(i)%mem_path_src, mname, pname)
      vm => get_virtual_model(mname)
      this%mvr(i)%is_provider_active = vm%is_local
      this%prov_is_m1(i) = associated(vm, this%model1)
      call split_mem_path(this%mvr(i)%mem_path_tgt, mname, pname)
      vm => get_virtual_model(mname)
      this%mvr(i)%is_receiver_active = vm%is_local
    end do

    ! loop over mvr's, if provider is active,
    ! store mapped feature index in array for sync
    allocate (pkg_mvr)

    do i = 1, nr_active_movers
      if (this%mvr(i)%is_provider_active) then
        ! store mapped feature id in array (for synchronization when parallel)
        call set_packagemover_pointer(pkg_mvr, this%mvr(i)%mem_path_src)
        if (this%prov_is_m1(i)) then
          this%id_mapped_m1(i) = pkg_mvr%iprmap(this%mvr(i)%iRchNrSrc)
          this%id_mapped_m2(i) = -1
        else
          this%id_mapped_m1(i) = -1
          this%id_mapped_m2(i) = pkg_mvr%iprmap(this%mvr(i)%iRchNrSrc)
        end if
      end if
    end do

  end subroutine xmvr_initialize_movers

  !> @brief Calculates qpactual and stores it for synchronization
  !<
  subroutine xmvr_cf(this)
    class(GwfExgMoverType) :: this
    ! local
    integer(I4B) :: i

    do i = 1, this%nmvr
      if (this%mvr(i)%is_provider_active) then

        call this%mvr(i)%update_provider()

        ! copy calculated rate to arrays for synchronization:
        if (this%prov_is_m1(i)) then
          this%qpactual_m1(i) = this%mvr(i)%qpactual
          this%qavailable_m1(i) = this%mvr(i)%qavailable
          this%qpactual_m2(i) = DNODATA
          this%qavailable_m2(i) = DNODATA
        else
          this%qpactual_m1(i) = DNODATA
          this%qavailable_m1(i) = DNODATA
          this%qpactual_m2(i) = this%mvr(i)%qpactual
          this%qavailable_m2(i) = this%mvr(i)%qavailable
        end if
      end if
    end do

  end subroutine xmvr_cf

  !> @brief Assign synced qpactual to mover and update receiver
  !<
  subroutine xmvr_fc(this)
    class(GwfExgMoverType) :: this
    ! local
    integer(I4B) :: i

    do i = 1, this%nmvr
      if (this%mvr(i)%is_receiver_active) then
        ! copy from synchronization arrays back into movers:
        if (this%prov_is_m1(i)) then
          this%mvr(i)%qpactual = this%qpactual_m1(i)
          this%mvr(i)%qavailable = this%qavailable_m1(i)
        else
          this%mvr(i)%qpactual = this%qpactual_m2(i)
          this%mvr(i)%qavailable = this%qavailable_m2(i)
        end if
        call this%mvr(i)%update_receiver()
      end if
    end do

  end subroutine xmvr_fc

  !> @brief Overrides budget routine to first assign the
  !< mapped features ids from the synchronization arrays
  subroutine xmvr_bd(this)
    class(GwfExgMoverType) :: this
    ! local
    integer(I4B) :: i

    ! copy from synchronization arrays back into movers:
    do i = 1, this%nmvr
      if (this%prov_is_m1(i)) then
        this%mvr(i)%iRchNrSrcMapped = this%id_mapped_m1(i)
      else
        this%mvr(i)%iRchNrSrcMapped = this%id_mapped_m2(i)
      end if
    end do

    call this%fill_budobj()

  end subroutine xmvr_bd

  subroutine xmvr_allocate_arrays(this)
    class(GwfExgMoverType) :: this
    ! local
    integer(I4B) :: i

    call this%GwfMvrType%allocate_arrays()

    allocate (this%prov_is_m1(this%maxmvr))
    call mem_allocate(this%qpactual_m1, this%maxmvr, 'QPACTUAL_M1', &
                      this%memoryPath)
    call mem_allocate(this%qpactual_m2, this%maxmvr, 'QPACTUAL_M2', &
                      this%memoryPath)
    call mem_allocate(this%qavailable_m1, this%maxmvr, 'QAVAILABLE_M1', &
                      this%memoryPath)
    call mem_allocate(this%qavailable_m2, this%maxmvr, 'QAVAILABLE_M2', &
                      this%memoryPath)
    call mem_allocate(this%id_mapped_m1, this%maxmvr, 'ID_MAPPED_M1', &
                      this%memoryPath)
    call mem_allocate(this%id_mapped_m2, this%maxmvr, 'ID_MAPPED_M2', &
                      this%memoryPath)

    do i = 1, this%maxmvr
      this%id_mapped_m1(i) = 0
      this%id_mapped_m2(i) = 0
      this%qpactual_m1(i) = DNODATA
      this%qpactual_m2(i) = DNODATA
      this%qavailable_m1(i) = DNODATA
      this%qavailable_m2(i) = DNODATA
    end do

  end subroutine xmvr_allocate_arrays

  subroutine xmvr_da(this)
    class(GwfExgMoverType) :: this

    call this%GwfMvrType%mvr_da()

    deallocate (this%prov_is_m1)
    call mem_deallocate(this%qpactual_m1)
    call mem_deallocate(this%qpactual_m2)
    call mem_deallocate(this%qavailable_m1)
    call mem_deallocate(this%qavailable_m2)
    call mem_deallocate(this%id_mapped_m1)
    call mem_deallocate(this%id_mapped_m2)

  end subroutine xmvr_da

end module GwfExgMoverModule
