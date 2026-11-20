!> @brief Time-varying hydraulic conductivity with grid array input (TVKG)
!!
!! Applies the K, K22 and K33 grid array values written by input in the
!! current time step to the NPF hydraulic conductivity arrays: values given
!! in a PERIOD block in rp, values written by a time-array series in ad.
!! Input arrays are indexed by user node; other values keep their last
!! setting.
!<
module TvkgModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENVARNAME, DZERO
  use BaseDisModule, only: DisBaseType
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error, count_errors, store_error_filename
  use SimVariablesModule, only: errmsg
  use TdisModule, only: kper, kstp
  use IdmParamMetaViewModule, only: IdmParamMetaViewType
  use IdmParamMetaModule, only: RS_VALUES, RS_SERIES, RS_NETCDF

  implicit none
  private
  public :: TvkgType
  public :: tvkg_cr

  !> @brief Time-varying hydraulic conductivity with grid array input
  !<
  type :: TvkgType
    character(len=LENMEMPATH) :: input_mempath = '' !< input context path
    character(len=LINELENGTH) :: input_fname = '' !< input file name
    character(len=LENMEMPATH) :: name_model = '' !< model name
    integer(I4B) :: iout = 0 !< listing file unit
    integer(I4B), pointer :: iprpak => null() !< print input flag
    class(DisBaseType), pointer :: dis => null() !< discretization
    integer(I4B), pointer :: ik22overk => null() !< NPF K22 is a ratio of K
    integer(I4B), pointer :: ik33overk => null() !< NPF K33 is a ratio of K
    real(DP), dimension(:), pointer, contiguous :: k11 => null() !< NPF K (before VSC adjustment if VSC is active)
    real(DP), dimension(:), pointer, contiguous :: k22 => null() !< NPF K22 (before VSC adjustment if VSC is active)
    real(DP), dimension(:), pointer, contiguous :: k33 => null() !< NPF K33 (before VSC adjustment if VSC is active)
    integer(I4B), pointer :: kchangeper => null() !< NPF period of last change
    integer(I4B), pointer :: kchangestp => null() !< NPF time step of last change
    integer(I4B), dimension(:), pointer, contiguous :: nodekchange => null() !< NPF node changed flags
    real(DP), dimension(:), pointer, contiguous :: k11_src => null() !< input K
    real(DP), dimension(:), pointer, contiguous :: k22_src => null() !< input K22
    real(DP), dimension(:), pointer, contiguous :: k33_src => null() !< input K33
    type(IdmParamMetaViewType) :: input_meta !< input loading information of K, K22 and K33
    integer(I4B), pointer :: iper => null() !< last period loaded
    integer(I4B) :: changed_kper = -1 !< stress period of the last change
    integer(I4B) :: changed_kstp = -1 !< time step of the last change
    logical(LGP) :: tas_active = .false. !< time-array series files given
  contains
    procedure :: ar
    procedure :: rp
    procedure :: ad
    procedure :: da
    procedure, private :: apply_written
    procedure, private :: check_netcdf_ratio
    procedure, private :: set_value
  end type TvkgType

contains

  !> @brief Create a TVKG object
  !<
  subroutine tvkg_cr(tvkg, name_model, mempath, input_fname, iout)
    type(TvkgType), pointer, intent(out) :: tvkg
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: mempath
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
    allocate (tvkg)
    tvkg%name_model = name_model
    tvkg%input_mempath = mempath
    tvkg%input_fname = input_fname
    tvkg%iout = iout
    allocate (tvkg%iprpak)
    tvkg%iprpak = 0
  end subroutine tvkg_cr

  !> @brief Set pointers to NPF and input arrays
  !<
  subroutine ar(this, dis)
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryManagerModule, only: get_isize
    class(TvkgType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    character(len=LENMEMPATH) :: npf_mempath
    logical(LGP) :: found
    integer(I4B), pointer :: invsc
    integer(I4B) :: isize
    character(len=*), parameter :: fmttvkg = &
      "(1x,/1x,'TVKG -- TIME-VARYING K GRID ARRAY INPUT, INPUT READ FROM &
      &MEMPATH ', a, //)"

    write (this%iout, fmttvkg) trim(this%input_mempath)
    this%dis => dis

    ! NPF arrays
    npf_mempath = create_mem_path(this%name_model, 'NPF')
    call mem_setptr(this%ik22overk, 'IK22OVERK', npf_mempath)
    call mem_setptr(this%ik33overk, 'IK33OVERK', npf_mempath)
    ! with VSC active, change the K values that VSC scales by the
    ! viscosity ratio each time step
    call mem_setptr(invsc, 'INVSC', npf_mempath)
    if (invsc > 0) then
      call mem_setptr(this%k11, 'K11INPUT', npf_mempath)
      call mem_setptr(this%k22, 'K22INPUT', npf_mempath)
      call mem_setptr(this%k33, 'K33INPUT', npf_mempath)
    else
      call mem_setptr(this%k11, 'K11', npf_mempath)
      call mem_setptr(this%k22, 'K22', npf_mempath)
      call mem_setptr(this%k33, 'K33', npf_mempath)
    end if
    call mem_setptr(this%kchangeper, 'KCHANGEPER', npf_mempath)
    call mem_setptr(this%kchangestp, 'KCHANGESTP', npf_mempath)
    call mem_setptr(this%nodekchange, 'NODEKCHANGE', npf_mempath)

    ! input arrays (user nodes)
    call mem_setptr(this%k11_src, 'K', this%input_mempath)
    call mem_setptr(this%k22_src, 'K22', this%input_mempath)
    call mem_setptr(this%k33_src, 'K33', this%input_mempath)

    ! options
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, found)
    call get_isize('TAS6_FILENAME', this%input_mempath, isize)
    this%tas_active = isize > 0

    ! input loading information
    call this%input_meta%init([character(len=LENVARNAME) :: 'K', 'K22', 'K33'], &
                              this%input_mempath)
    call mem_setptr(this%iper, 'IPER', this%input_mempath)
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine ar

  !> @brief Apply the values given in this period's PERIOD block
  !!
  !! Applies them once, at the time step the block is read.
  !<
  subroutine rp(this)
    class(TvkgType) :: this
    if (this%iper /= kper) return
    call this%check_netcdf_ratio()
    call this%apply_written(RS_VALUES)
  end subroutine rp

  !> @brief Apply the values written by time-array series in this time step
  !<
  subroutine ad(this)
    class(TvkgType) :: this
    if (.not. this%tas_active) return
    call this%apply_written(RS_SERIES)
  end subroutine ad

  !> @brief Apply the values written in this time step in the given way
  !!
  !! rs is RS_VALUES (given in a PERIOD block) or RS_SERIES (written by a
  !! time-array series). Change flags are reset at the first change in a
  !! time step, so they hold the changes from both rp and ad, and an ATS
  !! retry of the time step keeps them.
  !<
  subroutine apply_written(this, rs)
    class(TvkgType) :: this
    integer(I4B), intent(in) :: rs !< read state of the writes to apply
    logical(LGP), dimension(3) :: fields
    integer(I4B) :: j, nodeu, node

    call this%input_meta%update()
    do nodeu = 1, size(this%k11_src)
      do j = 1, 3
        fields(j) = this%input_meta%is_written(j, nodeu, rs)
      end do
      if (.not. any(fields)) cycle
      node = this%dis%get_nodenumber(nodeu, 0)
      if (node < 1) cycle
      if (this%changed_kper /= kper .or. this%changed_kstp /= kstp) then
        this%kchangeper = kper
        this%kchangestp = kstp
        this%nodekchange = 0
        this%changed_kper = kper
        this%changed_kstp = kstp
      end if
      ! K first, so K22 and K33 ratios use the new K
      if (fields(1)) &
        call this%set_value(node, nodeu, this%k11_src, this%k11, 'K', .false.)
      if (fields(2)) &
        call this%set_value(node, nodeu, this%k22_src, this%k22, 'K22', &
                            this%ik22overk == 1)
      if (fields(3)) &
        call this%set_value(node, nodeu, this%k33_src, this%k33, 'K33', &
                            this%ik33overk == 1)
    end do

    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine apply_written

  !> @brief Stop if a K22 or K33 ratio is read from NetCDF
  !!
  !! NetCDF period arrays hold the values in effect in each period, which
  !! a ratio of K does not.
  !<
  subroutine check_netcdf_ratio(this)
    class(TvkgType) :: this
    character(len=*), parameter :: fmterr = &
      "('TVKG ',a,' cannot be read from NetCDF input when it is a ratio of &
      &K (NPF ',a,'). Give ',a,' in the TVKG file.')"

    if (this%ik22overk == 1 .and. this%input_meta%read_state(2) == RS_NETCDF) then
      write (errmsg, fmterr) 'K22', 'K22OVERK', 'K22'
      call store_error(errmsg)
    end if
    if (this%ik33overk == 1 .and. this%input_meta%read_state(3) == RS_NETCDF) then
      write (errmsg, fmterr) 'K33', 'K33OVERK', 'K33'
      call store_error(errmsg)
    end if
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine check_netcdf_ratio

  !> @brief Set one NPF value from input
  !<
  subroutine set_value(this, node, nodeu, src, dst, name, ratio)
    class(TvkgType) :: this
    integer(I4B), intent(in) :: node
    integer(I4B), intent(in) :: nodeu
    real(DP), dimension(:), intent(in) :: src
    real(DP), dimension(:), intent(inout) :: dst
    character(len=*), intent(in) :: name
    logical(LGP), intent(in) :: ratio
    character(len=LINELENGTH) :: cellstr
    character(len=*), parameter :: fmtvalchg = &
      "(1x,'TVKG: setting ',a,' for cell ',a,' in period ',i0, &
      &' step ',i0,' = ',g12.5)"

    dst(node) = src(nodeu)
    if (ratio) dst(node) = dst(node) * this%k11(node)
    this%nodekchange(node) = 1

    if (dst(node) <= DZERO) then
      call this%dis%noder_to_string(node, cellstr)
      write (errmsg, '(a,1x,a,1x,a,1x,g15.6)') &
        'TVKG changed hydraulic property', name, &
        'is <= 0 for cell '//trim(cellstr), dst(node)
      call store_error(errmsg)
    end if
    if (this%iprpak /= 0) then
      call this%dis%noder_to_string(node, cellstr)
      write (this%iout, fmtvalchg) name, trim(cellstr), kper, kstp, dst(node)
    end if
  end subroutine set_value

  !> @brief Deallocate
  !<
  subroutine da(this)
    class(TvkgType) :: this
    deallocate (this%iprpak)
    nullify (this%dis, this%ik22overk, this%ik33overk)
    nullify (this%iper)
    nullify (this%k11, this%k22, this%k33)
    nullify (this%kchangeper, this%kchangestp, this%nodekchange)
    nullify (this%k11_src, this%k22_src, this%k33_src)
    call this%input_meta%destroy()
  end subroutine da

end module TvkgModule
