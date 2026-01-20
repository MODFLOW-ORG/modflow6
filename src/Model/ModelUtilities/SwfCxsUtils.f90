!> @brief This module contains stateless sfr subroutines and functions
!!
!! This module contains the functions to calculate the wetted perimeter
!! and cross-sectional area for a reach cross-section that are used in
!! the streamflow routing (SFR) package. It also contains subroutines to
!! calculate the wetted perimeter and cross-sectional area for each
!! line segment in the cross-section. This module does not depend on the
!! SFR package.
!!
!<
module SwfCxsUtilsModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DHALF, DTWOTHIRDS, DONE, DTWO, &
                             DEM5, DP999, DTHREE
  use SmoothingModule, only: sChSmooth

  implicit none
  private
  public :: calc_depth_from_q
  public :: calc_qman
  public :: get_saturated_topwidth
  public :: get_wetted_topwidth
  public :: get_wetted_perimeter
  public :: get_cross_section_area
  public :: get_hydraulic_radius
  public :: get_hydraulic_radius_xf
  public :: get_mannings_section
  public :: calc_composite_roughness
  public :: get_conveyance
  public :: get_composite_conveyance

contains

  !> @brief Calculate the depth at the midpoint of a irregular cross-section
  !<
  function calc_depth_from_q(qrch, width, rough, slope, &
                             cxs_xf, cxs_h, cxs_rf, unitconv, &
                             icalcmeth) result(depth)
    ! -- dummy variables
    real(DP), intent(in) :: qrch !< streamflow
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: unitconv !< conversion unit numerator in mannings equation
    integer(I4B), intent(in) :: icalcmeth !< manning calculation method

    ! -- return
    real(DP) :: depth !< stream depth at midpoint of reach

    if (icalcmeth == 1) then
      ! slower but robust bisection method
      depth = calc_depth_from_q_bisect(qrch, width, rough, slope, &
                                       cxs_xf, cxs_h, cxs_rf, unitconv, &
                                       icalcmeth)
    else
      ! faster but less forgiving newton-raphson method
      depth = calc_depth_from_q_nr(qrch, width, rough, slope, &
                                   cxs_xf, cxs_h, cxs_rf, unitconv, &
                                   icalcmeth)
    end if
  end function calc_depth_from_q

  !> @brief Calculate the depth at the midpoint of a irregular cross-section
  !!
  !! Method to calculate the depth at the midpoint of a reach with a
  !! irregular cross-section using bisection.
  !!
  !<
  function calc_depth_from_q_bisect(qrch, width, rough, slope, &
                                    cxs_xf, cxs_h, cxs_rf, unitconv, &
                                    icalcmeth) result(depth)
    ! -- dummy variables
    !class(SfrType) :: this !< SfrType object
    !integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: qrch !< streamflow
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: unitconv !< conversion unit numerator in mannings equation
    integer(I4B), intent(in) :: icalcmeth !< manning calculation method

    ! -- return
    real(DP) :: depth !< stream depth at midpoint of reach

    ! -- local variables
    integer(I4B) :: maxiter = 100
    integer(I4B) :: iter
    real(DP) :: dmaxchg = DEM5
    real(DP) :: a
    real(DP) :: b
    real(DP) :: c
    real(DP) :: f_a
    real(DP) :: f_c

    ! constrain the bisection range
    a = DZERO
    b = maxval(cxs_h) * DTWO
    !
    ! -- bisection iteration
    bisectiter: do iter = 1, maxiter

      c = (a + b) / DTWO
      f_c = calc_qman(c, width, rough, slope, &
                      cxs_xf, cxs_h, cxs_rf, unitconv, icalcmeth) - qrch
      if (f_c == DZERO .or. (b - a) / DTWO < dmaxchg) then
        depth = c
        return
      end if

      f_a = calc_qman(a, width, rough, slope, &
                      cxs_xf, cxs_h, cxs_rf, unitconv, icalcmeth) - qrch
      if (sign(DONE, f_c) == sign(DONE, f_a)) then
        a = c
      else
        b = c
      end if

    end do bisectiter
    !
    ! TODO: raise an error
    print *, "bisection routine failed"
  end function calc_depth_from_q_bisect

  !> @brief Calculate the depth at the midpoint of a irregular cross-section
  !!
  !! Method to calculate the depth at the midpoint of a reach with a
  !! irregular cross-section using Newton-Raphson.
  !!
  !<
  function calc_depth_from_q_nr(qrch, width, rough, slope, &
                                cxs_xf, cxs_h, cxs_rf, unitconv, &
                                icalcmeth) result(depth)
    ! -- dummy variables
    !class(SfrType) :: this !< SfrType object
    !integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: qrch !< streamflow
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: unitconv !< conversion unit numerator in mannings equation
    integer(I4B), intent(in) :: icalcmeth !< manning calculation method

    ! -- return
    real(DP) :: depth !< stream depth at midpoint of reach

    ! -- local variables
    integer(I4B) :: maxiter = 100
    integer(I4B) :: iter
    real(DP) :: dmaxchg = DEM5
    real(DP) :: deps = DEM5 * DP999
    real(DP) :: perturbation
    real(DP) :: q0
    real(DP) :: q1
    real(DP) :: dq
    real(DP) :: derv
    real(DP) :: d
    real(DP) :: dd
    real(DP) :: residual
    !
    ! -- initialize variables
    perturbation = deps * DTWO
    d = DZERO
    q0 = DZERO
    residual = q0 - qrch
    !
    ! -- Newton-Raphson iteration
    nriter: do iter = 1, maxiter
      !call this%sfr_calc_qman(n, d + perturbation, q1)
      q1 = calc_qman(d + perturbation, width, rough, slope, &
                     cxs_xf, cxs_h, cxs_rf, unitconv, icalcmeth)
      dq = (q1 - q0)
      if (dq /= DZERO) then
        derv = perturbation / (q1 - q0)
      else
        derv = DZERO
      end if
      dd = derv * residual
      d = d - dd
      !call this%sfr_calc_qman(n, d, q0)
      q0 = calc_qman(d, width, rough, slope, &
                     cxs_xf, cxs_h, cxs_rf, unitconv, icalcmeth)
      residual = q0 - qrch
      !
      ! -- check for convergence
      if (abs(dd) < dmaxchg) then
        exit nriter
      end if

    end do nriter
    depth = d
  end function calc_depth_from_q_nr

  !> @brief Calculate streamflow using Manning's equation
  !!
  !! Method to calculate the streamflow using Manning's equation for a
  !! single reach defined by a cross section.
  !!
  !<
  function calc_qman(depth, width, rough, slope, &
                     cxs_xf, cxs_h, cxs_rf, unitconv, icalcmeth) result(qman)

    ! -- dummy variables
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: unitconv !< conversion unit numerator in mannings equation
    integer(I4B), intent(in) :: icalcmeth ! 0 is composite linear, 1 is by section, 2 is composite nonlinear

    ! return value
    real(DP) :: qman !< calculated mannings flow

    ! -- local variables
    integer(I4B) :: linmeth

    select case (icalcmeth)
    case (0) ! composite linear
      linmeth = 0
      qman = calc_qman_composite(depth, width, rough, slope, &
                                 cxs_xf, cxs_h, cxs_rf, unitconv, &
                                 linmeth)
    case (1) ! by section
      qman = calc_qman_by_section(depth, width, rough, slope, &
                                  cxs_xf, cxs_h, cxs_rf, unitconv)
    case (2) ! composite nonlinear
      linmeth = 1
      qman = calc_qman_composite(depth, width, rough, slope, &
                                 cxs_xf, cxs_h, cxs_rf, unitconv, &
                                 linmeth)
    end select
  end function calc_qman

  function calc_qman_composite(depth, width, rough, slope, &
                               cxs_xf, cxs_h, cxs_rf, unitconv, &
                               linmeth) result(qman)

    ! -- dummy variables
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: unitconv !< conversion unit numerator in mannings equation
    integer(I4B), intent(in) :: linmeth !< method for composite (0) linear (1) nonlinear

    ! return value
    real(DP) :: qman !< calculated mannings flow

    ! -- local variables
    integer(I4B) :: npts
    real(DP) :: sat
    real(DP) :: derv
    real(DP) :: aw
    real(DP) :: wp
    real(DP) :: rh
    real(DP) :: rough_composite
    !
    ! -- initialize variables
    qman = DZERO
    !
    ! -- calculate Manning's discharge for non-zero depths
    if (depth > DZERO) then

      npts = size(cxs_xf)
      rough_composite = calc_composite_roughness(npts, depth, width, &
                                                 rough, slope, &
                                                 cxs_xf, cxs_h, cxs_rf, &
                                                 linmeth)
      wp = get_wetted_perimeter(npts, cxs_xf, cxs_h, width, depth)
      aw = get_cross_section_area(npts, cxs_xf, cxs_h, width, depth)
      if (wp > DZERO) then
        rh = aw / wp
      else
        rh = DZERO
      end if
      qman = unitconv * aw * (rh**DTWOTHIRDS) * sqrt(slope) / rough
      call sChSmooth(depth, sat, derv)
      qman = sat * qman
    end if
  end function calc_qman_composite

  function calc_composite_roughness(npts, depth, width, rough, slope, &
                                    cxs_xf, cxs_h, cxs_rf, &
                                    linmeth) result(rc)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    integer(I4B), intent(in) :: linmeth !< method for composite calculation; linear 0 or nonlinear 1

    ! return value
    real(DP) :: rc !< composite roughness

    ! -- local variables
    integer(I4B) :: n
    real(DP) :: wp
    real(DP) :: exp1
    real(DP) :: exp2
    real(DP), dimension(npts) :: stations
    real(DP), dimension(npts - 1) :: perimeters
    !
    ! stations
    do n = 1, npts
      stations(n) = cxs_xf(n) * width
    end do
    !
    ! -- select method
    select case (linmeth)
    case (0) ! linear
      exp1 = DONE
      exp2 = DONE
    case (1) ! nonlinear
      exp1 = 1.5D0
      exp2 = DTWO / DTHREE
    end select
    !
    ! -- initialize variables
    rc = rough
    !
    ! -- calculate composite roughness
    if (depth > DZERO) then

      if (npts > 1) then
        call get_wetted_perimeters(npts, stations, cxs_h, depth, perimeters)
        wp = DZERO
        rc = DZERO
        do n = 1, npts - 1
          rc = rc + (rough * cxs_rf(n))**exp1 * perimeters(n)
          wp = wp + perimeters(n)
        end do
        if (wp > DZERO) then
          rc = (rc / wp)**exp2
        else
          rc = rough
        end if
      end if

    end if
  end function calc_composite_roughness

  function calc_qman_by_section(depth, width, rough, slope, &
                                cxs_xf, cxs_h, cxs_rf, unitconv) result(qman)

    ! -- dummy variables
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< mannings reach roughness coefficient
    real(DP), intent(in) :: slope !< reach bottom slope
    real(DP), dimension(:), intent(in) :: cxs_xf ! xfraction distances for this cross section
    real(DP), dimension(:), intent(in) :: cxs_h ! heights for this cross section
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: unitconv !< conversion unit numerator in mannings equation

    ! return value
    real(DP) :: qman !< calculated mannings flow

    ! -- local variables
    integer(I4B) :: npts
    real(DP) :: sat
    real(DP) :: derv
    real(DP) :: aw
    real(DP) :: wp
    real(DP) :: rh
    !
    ! -- initialize variables
    qman = DZERO
    !
    ! -- calculate Manning's discharge for non-zero depths
    if (depth > DZERO) then
      npts = size(cxs_xf)
      !
      ! -- set constant terms for Manning's equation
      call sChSmooth(depth, sat, derv)
      !
      ! -- calculate the mannings coefficient that is a
      !    function of depth
      if (npts > 1) then
        !
        ! -- call function that calculates flow for an
        !    n-point cross section
        qman = get_mannings_section(npts, &
                                    cxs_xf, &
                                    cxs_h, &
                                    cxs_rf, &
                                    rough, &
                                    unitconv, &
                                    width, &
                                    slope, &
                                    depth)
      else

        ! hydraulically wide channel (only 1 point is defined)
        aw = width * depth
        wp = width
        if (wp > DZERO) then
          rh = aw / wp
        else
          rh = DZERO
        end if
        qman = unitconv * aw * (rh**DTWOTHIRDS) * sqrt(slope) / rough
      end if
      !
      ! -- calculate stream flow
      qman = sat * qman
    end if
  end function calc_qman_by_section

  !> @brief Calculate the saturated top width for a reach
  !!
  !! Function to calculate the maximum top width for a reach using the
  !! cross-section station data.
  !!
  !! @return      w               saturated top width
  !<
  function get_saturated_topwidth(npts, xfraction, width) result(w)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station fractional distances (x-distance)
    real(DP), intent(in) :: width
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: w
    real(DP), dimension(npts) :: stations
    !
    ! -- calculate station from xfractions and width
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- calculate the saturated top width
    if (npts > 1) then
      w = stations(npts) - stations(1)
    else
      w = stations(1)
    end if
  end function get_saturated_topwidth

  !> @brief Calculate the wetted top width for a reach
  !!
  !! Function to calculate the wetted top width for a reach using the
  !! cross-section station-height data given a passed depth.
  !!
  !! @return      w               wetted top width
  !<
  function get_wetted_topwidth(npts, xfraction, heights, width, d) result(w)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: w
    real(DP), dimension(npts) :: stations
    real(DP), dimension(npts - 1) :: widths
    !
    ! -- calculate station from xfractions and width
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- initialize the wetted perimeter for the reach
    w = DZERO
    !
    ! -- calculate the wetted top width for each line segment
    call get_wetted_topwidths(npts, stations, heights, d, widths)
    !
    ! -- calculate the wetted top widths
    do n = 1, npts - 1
      w = w + widths(n)
    end do
  end function get_wetted_topwidth

  !> @brief Calculate the wetted perimeter for a reach
  !!
  !! Function to calculate the wetted perimeter for a reach using the
  !! cross-section station-height data given a passed depth.
  !!
  !! @return      p               wetted perimeter
  !<
  function get_wetted_perimeter(npts, xfraction, heights, width, d) result(p)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section x fractions
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: width !< cross section width
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: p
    real(DP), dimension(npts) :: stations
    real(DP), dimension(npts - 1) :: perimeters
    !
    ! -- initialize stations
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- initialize the wetted perimeter for the reach
    p = DZERO
    !
    ! -- calculate the wetted perimeter for each line segment
    call get_wetted_perimeters(npts, stations, heights, d, perimeters)
    !
    ! -- calculate the wetted perimenter
    do n = 1, npts - 1
      p = p + perimeters(n)
    end do
  end function get_wetted_perimeter

  !> @brief Calculate the cross-sectional area for a reach
  !!
  !! Function to calculate the cross-sectional area for a reach using
  !! the cross-section station-height data given a passed depth.
  !!
  !! @return      a               cross-sectional area
  !<
  function get_cross_section_area(npts, xfraction, heights, width, d) result(a)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: a
    real(DP), dimension(npts) :: stations
    real(DP), dimension(npts - 1) :: areas
    !
    ! -- calculate station from xfractions and width
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- initialize the area
    a = DZERO
    !
    ! -- calculate the cross-sectional area for each line segment
    call get_cross_section_areas(npts, stations, heights, d, areas)
    !
    ! -- calculate the cross-sectional area
    do n = 1, npts - 1
      a = a + areas(n)
    end do
  end function get_cross_section_area

  !> @brief Calculate conveyance
  !!
  !! Return a calculated conveyance.  Calculation
  !! depends on whether the channel is a rectangular channel
  !< with vertical sides or something else.
  function get_conveyance(npts, xfraction, heights, cxs_rf, &
                          width, rough, d) result(c)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< reach roughness
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! result
    real(DP) :: c
    if (is_rectangular(xfraction)) then
      c = get_rectangular_conveyance(npts, xfraction, heights, cxs_rf, &
                                     width, rough, d)
    else
      c = get_composite_conveyance(npts, xfraction, heights, cxs_rf, &
                                   width, rough, d)
    end if
  end function get_conveyance

  !> @brief Calculate conveyance for rectangular channel
  !<
  function get_rectangular_conveyance(npts, xfraction, heights, cxs_rf, &
                                      width, rough, d) result(c)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< reach roughness
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    real(DP) :: c
    real(DP) :: a
    real(DP) :: p
    real(DP) :: ravg

    ! find cross sectional area and wetted perimeter
    a = get_cross_section_area(npts, xfraction, heights, width, d)
    p = get_wetted_perimeter(npts, xfraction, heights, width, d)

    ! calculate roughness or determine an average
    if (has_uniform_resistance(cxs_rf)) then
      ravg = rough * cxs_rf(1)
    else
      ravg = calc_composite_roughness(npts, d, width, rough, DZERO, &
                                      xfraction, heights, cxs_rf, 0)
    end if

    ! make conveyance calculation
    c = DZERO
    if (p > DZERO) then
      c = a / ravg * (a / p)**DTWOTHIRDS
    end if
  end function get_rectangular_conveyance

  !> @brief Determine if cross section is rectangular
  !<
  function is_rectangular(xfraction)
    ! dummy
    real(DP), dimension(:), intent(in) :: xfraction !< cross-section station distances (x-distance)
    ! result
    logical(LGP) :: is_rectangular
    is_rectangular = .false.
    if (size(xfraction) == 4) then
      if (xfraction(1) == xfraction(2) .and. &
          xfraction(3) == xfraction(4)) then
        is_rectangular = .true.
      end if
    end if
  end function is_rectangular

  !> @brief Determine if roughness is uniform for the section
  !<
  function has_uniform_resistance(cxs_rf)
    ! dummy
    real(DP), dimension(:), intent(in) :: cxs_rf !< cross-section station distances (x-distance)
    ! result
    logical(LGP) :: has_uniform_resistance
    ! local
    real(DP) :: rmin
    real(DP) :: rmax
    rmin = minval(cxs_rf(1:3))
    rmax = maxval(cxs_rf(1:3))
    has_uniform_resistance = (rmin == rmax)
  end function has_uniform_resistance

  !> @brief Calculate composite conveyance
  !!
  !! Function to calculate the composite conveyance.
  !! This is based on the approach in HEC-RAS, where
  !! a conveyance is calculated for each line segment
  !! in the cross section, and then summed to produce
  !! a total conveyance.
  !!
  !<
  function get_composite_conveyance(npts, xfraction, heights, cxs_rf, &
                                    width, rough, d) result(c)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), dimension(:), intent(in) :: cxs_rf ! mannings fractions for this cross section
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: rough !< reach roughness
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: c
    real(DP) :: p
    real(DP) :: rc
    real(DP) :: rh
    real(DP) :: cn
    real(DP), dimension(npts) :: stations
    real(DP), dimension(npts - 1) :: areas
    real(DP), dimension(npts - 1) :: perimeters
    !
    ! -- calculate station from xfractions and width
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- calculate the cross-sectional area for each line segment
    call get_cross_section_areas(npts, stations, heights, d, areas)
    !
    ! -- calculate the wetted perimeter for each line segment
    call get_wetted_perimeters(npts, stations, heights, d, perimeters)
    !
    ! -- calculate the composite conveyance
    c = DZERO
    do n = 1, npts - 1
      rc = cxs_rf(n) * rough
      p = perimeters(n)
      cn = DZERO
      if (p > DZERO) then
        rh = areas(n) / p
        cn = areas(n) / rc * rh**DTWOTHIRDS
      end if
      c = c + cn
    end do
  end function get_composite_conveyance

  !> @brief Calculate the hydraulic radius for a reach
  !!
  !! Function to calculate the hydraulic radius for a reach using
  !! the cross-section station-height data given a passed depth.
  !!
  !! @return      r               hydraulic radius
  !<
  function get_hydraulic_radius(npts, stations, heights, d) result(r)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: stations !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: r
    real(DP) :: p
    real(DP) :: a
    real(DP), dimension(npts - 1) :: areas
    real(DP), dimension(npts - 1) :: perimeters
    !
    ! -- initialize the hydraulic radius, perimeter, and area
    r = DZERO
    p = DZERO
    a = DZERO
    !
    ! -- calculate the wetted perimeter for each line segment
    call get_wetted_perimeters(npts, stations, heights, d, perimeters)
    !
    ! -- calculate the wetted perimenter
    do n = 1, npts - 1
      p = p + perimeters(n)
    end do
    !
    ! -- calculate the hydraulic radius only if the perimeter is non-zero
    if (p > DZERO) then
      !
      ! -- calculate the cross-sectional area for each line segment
      call get_cross_section_areas(npts, stations, heights, d, areas)
      !
      ! -- calculate the cross-sectional area
      do n = 1, npts - 1
        a = a + areas(n)
      end do
      !
      ! -- calculate the hydraulic radius
      r = a / p
    end if
  end function get_hydraulic_radius

  !> @brief Calculate the hydraulic radius for a reach
  !!
  !! Function to calculate the hydraulic radius for a reach using
  !! the cross-section xfraction-height data given a passed depth.
  !! This is different from get_hydraulic_radius as it requires
  !! xfraction and width instead of station.
  !!
  !! @return      r               hydraulic radius
  !<
  function get_hydraulic_radius_xf(npts, xfraction, heights, width, d) result(r)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: r
    real(DP), dimension(npts) :: stations
    !
    ! -- calculate station from xfractions and width
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- calculate hydraulic radius
    r = get_hydraulic_radius(npts, stations, heights, d)
  end function get_hydraulic_radius_xf

  !> @brief Calculate the manning's discharge for a reach
  !!
  !! Function to calculate the mannings discharge for a reach
  !! by calculating the discharge for each section, which can
  !! have a unique Manning's coefficient given a passed depth.
  !!
  !! @return      q               reach discharge
  !<
  function get_mannings_section(npts, xfraction, heights, roughfracs, &
                                roughness, conv_fact, width, slope, d) result(q)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: xfraction !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), dimension(npts), intent(in) :: roughfracs !< cross-section Mannings roughness fraction data
    real(DP), intent(in) :: roughness !< base reach roughness
    real(DP), intent(in) :: conv_fact !< unit conversion factor
    real(DP), intent(in) :: width !< reach width
    real(DP), intent(in) :: slope !< reach slope
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: q
    real(DP) :: rh
    real(DP) :: r
    real(DP) :: p
    real(DP) :: a
    real(DP), dimension(npts) :: stations
    real(DP), dimension(npts - 1) :: areas
    real(DP), dimension(npts - 1) :: perimeters
    !
    ! -- initialize the hydraulic radius, perimeter, and area
    q = DZERO
    rh = DZERO
    r = DZERO
    p = DZERO
    a = DZERO
    !
    ! -- calculate station from xfractions and width
    do n = 1, npts
      stations(n) = xfraction(n) * width
    end do
    !
    ! -- calculate the wetted perimeter for each line segment
    call get_wetted_perimeters(npts, stations, heights, d, perimeters)
    !
    ! -- calculate the wetted perimenter
    do n = 1, npts - 1
      p = p + perimeters(n)
    end do
    !
    ! -- calculate the hydraulic radius only if the perimeter is non-zero
    if (p > DZERO) then
      !
      ! -- calculate the cross-sectional area for each line segment
      call get_cross_section_areas(npts, stations, heights, d, areas)
      !
      ! -- calculate the cross-sectional area
      do n = 1, npts - 1
        p = perimeters(n)
        r = roughness * roughfracs(n)
        if (p * r > DZERO) then
          a = areas(n)
          rh = a / p
          q = q + conv_fact * a * rh**DTWOTHIRDS * sqrt(slope) / r
        end if
      end do
    end if
  end function get_mannings_section

  ! -- private functions and subroutines

  !> @brief Calculate the wetted perimeters for each line segment
  !!
  !! Subroutine to calculate the wetted perimeter for each line segment
  !! that defines the reach using the cross-section station-height
  !! data given a passed depth.
  !!
  !<
  subroutine get_wetted_perimeters(npts, stations, heights, d, p)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: stations !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    real(DP), dimension(npts - 1), intent(inout) :: p !< wetted perimeter for each line segment
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: x0
    real(DP) :: x1
    real(DP) :: d0
    real(DP) :: d1
    real(DP) :: dmax
    real(DP) :: dmin
    real(DP) :: xlen
    real(DP) :: dlen
    !
    ! -- iterate over the station-height data
    do n = 1, npts - 1
      !
      ! -- initialize the wetted perimeter
      p(n) = DZERO
      !
      ! -- initialize station-height data for segment
      x0 = stations(n)
      x1 = stations(n + 1)
      d0 = heights(n)
      d1 = heights(n + 1)
      !
      ! -- get the start and end station position of the wetted segment
      call get_wetted_station(x0, x1, d0, d1, dmax, dmin, d)
      !
      ! -- calculate the wetted perimeter for the segment
      xlen = x1 - x0
      dlen = DZERO
      if (xlen > DZERO) then
        if (d > dmax) then
          dlen = dmax - dmin
        else
          dlen = d - dmin
        end if
      else
        if (d > dmin) then
          dlen = min(d, dmax) - dmin
        else
          dlen = DZERO
        end if
      end if
      p(n) = sqrt(xlen**DTWO + dlen**DTWO)
    end do
  end subroutine get_wetted_perimeters

  !> @brief Calculate the cross-sectional areas for each line segment
  !!
  !! Subroutine to calculate the cross-sectional area for each line segment
  !! that defines the reach using the cross-section station-height
  !! data given a passed depth.
  !!
  !<
  subroutine get_cross_section_areas(npts, stations, heights, d, a)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: stations !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    real(DP), dimension(npts - 1), intent(inout) :: a !< cross-sectional area for each line segment
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: x0
    real(DP) :: x1
    real(DP) :: d0
    real(DP) :: d1
    real(DP) :: dmax
    real(DP) :: dmin
    real(DP) :: xlen
    !
    ! -- iterate over the station-height data
    do n = 1, npts - 1
      !
      ! -- initialize the cross-sectional area
      a(n) = DZERO
      !
      ! -- initialize station-height data for segment
      x0 = stations(n)
      x1 = stations(n + 1)
      d0 = heights(n)
      d1 = heights(n + 1)
      !
      ! -- get the start and end station position of the wetted segment
      call get_wetted_station(x0, x1, d0, d1, dmax, dmin, d)
      !
      ! -- calculate the cross-sectional area for the segment
      xlen = x1 - x0
      if (xlen > DZERO) then
        !
        ! -- add the area above dmax
        if (d > dmax) then
          a(n) = xlen * (d - dmax)
        end if
        !
        ! -- add the area below dmax
        if (dmax /= dmin .and. d > dmin) then
          if (d < dmax) then
            a(n) = a(n) + DHALF * (d - dmin) * xlen
          else
            a(n) = a(n) + DHALF * (dmax - dmin) * xlen
          end if
        end if
      end if
    end do
  end subroutine get_cross_section_areas

  !> @brief Calculate the wetted top widths for each line segment
  !!
  !! Subroutine to calculate the wetted top width for each line segment
  !! that defines the reach using the cross-section station-height
  !! data given a passed depth.
  !!
  !<
  subroutine get_wetted_topwidths(npts, stations, heights, d, w)
    ! -- dummy variables
    integer(I4B), intent(in) :: npts !< number of station-height data for a reach
    real(DP), dimension(npts), intent(in) :: stations !< cross-section station distances (x-distance)
    real(DP), dimension(npts), intent(in) :: heights !< cross-section height data
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    real(DP), dimension(npts - 1), intent(inout) :: w !< wetted top widths for each line segment
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: x0
    real(DP) :: x1
    real(DP) :: d0
    real(DP) :: d1
    real(DP) :: dmax
    real(DP) :: dmin
    !
    ! -- iterate over the station-height data
    do n = 1, npts - 1
      !
      ! -- initialize station-height data for segment
      x0 = stations(n)
      x1 = stations(n + 1)
      d0 = heights(n)
      d1 = heights(n + 1)
      !
      ! -- get the start and end station position of the wetted segment
      call get_wetted_station(x0, x1, d0, d1, dmax, dmin, d)
      !
      ! -- calculate the wetted top width for the segment
      w(n) = x1 - x0
    end do
  end subroutine get_wetted_topwidths

  !> @brief Calculate the station values for the wetted portion of the cross-section
  !!
  !! Subroutine to calculate the station values that define the extent of the
  !! wetted portion of the cross section for a line segment. The left (x0) and
  !! right (x1) station positions are altered if the passed depth is less
  !! than the maximum line segment depth. If the line segment is dry the left
  !! and right station are equal. Otherwise the wetted station values are equal
  !! to the full line segment or smaller if the passed depth is less than
  !! the maximum line segment depth.
  !!
  !<
  pure subroutine get_wetted_station(x0, x1, d0, d1, dmax, dmin, d)
    ! -- dummy variables
    real(DP), intent(inout) :: x0 !< left station position
    real(DP), intent(inout) :: x1 !< right station position
    real(DP), intent(in) :: d0 !< depth at the left station
    real(DP), intent(in) :: d1 !< depth at the right station
    real(DP), intent(inout) :: dmax !< maximum depth
    real(DP), intent(inout) :: dmin !< minimum depth
    real(DP), intent(in) :: d !< depth to evaluate cross-section
    ! -- local variables
    real(DP) :: xlen
    real(DP) :: dlen
    real(DP) :: slope
    real(DP) :: xt
    !
    ! -- calculate the minimum and maximum depth
    dmin = min(d0, d1)
    dmax = max(d0, d1)
    !
    ! -- calculate x0 and x1
    if (d <= dmin) then
      ! -- if d is less than or equal to the minimum value the
      !    station length (xlen) is zero
      x1 = x0
    else if (d >= dmax) then
      ! x0 and x1 unchanged (full wetted width)
      continue
    else
      ! -- if d is between dmin and dmax, station length is less
      !    than d1 - d0
      xlen = x1 - x0
      dlen = d1 - d0
      ! because of preceding checks
      ! we know dmin<d<dmax, dmax>dmin, dlen > 0
      ! no need to check for dlen == 0
      slope = xlen / dlen
      xt = x0 + slope * (d - d0)
      if (d0 > d1) then
        ! x1 unchanged
        x0 = xt
      else
        !x0 unchanged
        x1 = xt
      end if
    end if
  end subroutine get_wetted_station

end module SwfCxsUtilsModule
