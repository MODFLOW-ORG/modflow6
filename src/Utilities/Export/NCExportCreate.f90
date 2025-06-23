!> @brief This module contains the NCExportCreateModule
!!
!! This module creates derived model netcdf export
!! objects. It is dependent on netcdf libraries.
!!
!<
module NCExportCreateModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: DIS, DISU, DISV
  use SimModule, only: store_error, store_error_filename
  use NumericalModelModule, only: NumericalModelType
  use BaseDisModule, only: DisBaseType
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use DisuModule, only: DisuType
  use ModelExportModule, only: export_models, get_export_model
  use ModelExportModule, only: ExportModelType
  use NCModelExportModule, only: ExportPackageType

  implicit none
  private
  public :: nc_export_create

contains

  !> @brief create model netcdf export type
  !!
  subroutine create_nc_export(export_model, num_model)
    use NCModelExportModule, only: NETCDF_MESH2D, NETCDF_STRUCTURED
    use MeshDisModelModule, only: Mesh2dDisExportType
    use MeshDisvModelModule, only: Mesh2dDisvExportType
    use DisNCStructuredModule, only: DisNCStructuredType
    use InputLoadTypeModule, only: ModelDynamicPkgsType
    type(ExportModelType), pointer, intent(inout) :: export_model
    class(NumericalModelType), pointer, intent(in) :: num_model
    class(Mesh2dDisExportType), pointer :: ugrid_dis
    class(Mesh2dDisvExportType), pointer :: ugrid_disv
    class(DisNCStructuredType), pointer :: structured_dis
    class(DisBaseType), pointer :: disbase

    select case (export_model%disenum)
    case (DIS)
      ! allocate nc structured grid export object
      if (export_model%nctype == NETCDF_MESH2D) then
        ! allocate nc structured grid export object
        allocate (ugrid_dis)

        ! set dis base type
        disbase => num_model%dis
        select type (disbase)
        type is (DisType)
          ugrid_dis%dis => disbase
        end select

        ! set dynamic loaders
        call create_export_pkglist(ugrid_dis%pkglist, export_model%loaders, &
                                   export_model%iout)

        ! initialize export object
        call ugrid_dis%init(export_model%modelname, export_model%modeltype, &
                            export_model%modelfname, export_model%nc_fname, &
                            export_model%disenum, NETCDF_MESH2D, &
                            export_model%iout)

        ! define export object
        call ugrid_dis%df()

        ! set base pointer
        export_model%nc_export => ugrid_dis
      else if (export_model%nctype == NETCDF_STRUCTURED) then
        ! allocate nc structured grid export object
        allocate (structured_dis)

        ! set dis base type
        disbase => num_model%dis
        select type (disbase)
        type is (DisType)
          structured_dis%dis => disbase
        end select

        ! set dynamic loaders
        call create_export_pkglist(structured_dis%pkglist, export_model%loaders, &
                                   export_model%iout)

        ! initialize export object
        call structured_dis%init(export_model%modelname, export_model%modeltype, &
                                 export_model%modelfname, export_model%nc_fname, &
                                 export_model%disenum, NETCDF_STRUCTURED, &
                                 export_model%iout)

        ! define export object
        call structured_dis%df()

        ! set base pointer
        export_model%nc_export => structured_dis
      end if
    case (DISV)
      if (export_model%nctype == NETCDF_MESH2D) then
        ! allocate nc structured grid export object
        allocate (ugrid_disv)

        ! set dis base type
        disbase => num_model%dis
        select type (disbase)
        type is (DisvType)
          ugrid_disv%disv => disbase
        end select

        ! set dynamic loaders
        call create_export_pkglist(ugrid_disv%pkglist, export_model%loaders, &
                                   export_model%iout)

        ! initialize export object
        call ugrid_disv%init(export_model%modelname, export_model%modeltype, &
                             export_model%modelfname, export_model%nc_fname, &
                             export_model%disenum, NETCDF_MESH2D, &
                             export_model%iout)

        ! define export object
        call ugrid_disv%df()

        ! set base pointer
        export_model%nc_export => ugrid_disv
      else
        errmsg = 'DISV model discretization only &
                 &supported as UGRID NetCDF export. &
                 &Model='//trim(export_model%modelname)//'.'
        call store_error(errmsg)
        call store_error_filename(export_model%modelfname)
      end if
    case default
      errmsg = 'Unsupported discretization for NetCDF model export. &
               &Model='//trim(export_model%modelname)//'.'
      call store_error(errmsg)
      call store_error_filename(export_model%modelfname)
    end select
  end subroutine create_nc_export

  subroutine create_export_pkglist(pkglist, loaders, iout)
    use ListModule, only: ListType
    use MemoryManagerExtModule, only: mem_set_value
    use InputLoadTypeModule, only: ModelDynamicPkgsType
    use InputLoadTypeModule, only: DynamicPkgLoadBaseType
    use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
    use LayerArrayLoadModule, only: LayerArrayLoadType
    use GridArrayLoadModule, only: GridArrayLoadType
    use IdmMf6FileModule, only: Mf6FileDynamicPkgLoadType
    type(ListType), intent(inout) :: pkglist
    type(ModelDynamicPkgsType), pointer, intent(in) :: loaders
    integer(I4B), intent(in) :: iout
    class(DynamicPkgLoadBaseType), pointer :: dynamic_pkg
    class(AsciiDynamicPkgLoadBaseType), pointer :: rp_loader
    type(ExportPackageType), pointer :: export_pkg
    integer(I4B), pointer :: export_arrays
    class(*), pointer :: obj
    logical(LGP) :: found, readasarrays
    integer(I4B) :: n

    ! create list of in scope loaders
    allocate (export_arrays)

    do n = 1, loaders%pkglist%Count()
      ! initialize export arrays option
      export_arrays = 0

      dynamic_pkg => loaders%get(n)

      ! update export arrays option
      call mem_set_value(export_arrays, 'EXPORT_NC', &
                         dynamic_pkg%mf6_input%mempath, found)

      readasarrays = (dynamic_pkg%readarraylayer .or. dynamic_pkg%readarraygrid)
      if (export_arrays > 0 .and. readasarrays) then
        select type (dynamic_pkg)
        type is (Mf6FileDynamicPkgLoadType)
          rp_loader => dynamic_pkg%rp_loader
          select type (rp_loader)
          type is (LayerArrayLoadType)
            ! create the export object
            allocate (export_pkg)
            call export_pkg%init(rp_loader%mf6_input, &
                                 rp_loader%bound_context%mshape, &
                                 rp_loader%bound_context%naux, &
                                 rp_loader%param_names, rp_loader%nparam)
            obj => export_pkg
            call pkglist%add(obj)
          type is (GridArrayLoadType)
            ! create the export object
            allocate (export_pkg)
            call export_pkg%init(rp_loader%mf6_input, &
                                 rp_loader%bound_context%mshape, &
                                 rp_loader%bound_context%naux, &
                                 rp_loader%param_names, rp_loader%nparam)
            obj => export_pkg
            call pkglist%add(obj)
          end select
        end select
      end if
    end do

    ! cleanup
    deallocate (export_arrays)
  end subroutine create_export_pkglist

  !> @brief initialize netcdf model export type
  !!
  subroutine nc_export_create()
    use NumericalModelModule, only: GetNumericalModelFromList
    use ListsModule, only: basemodellist
    use NCModelExportModule, only: NETCDF_UNDEF
    integer(I4B) :: n
    type(ExportModelType), pointer :: export_model
    class(NumericalModelType), pointer :: num_model
    integer(I4B) :: im
    do n = 1, export_models%Count()
      ! set pointer to export model
      export_model => get_export_model(n)
      if (export_model%nctype /= NETCDF_UNDEF) then
        ! netcdf export is active identify model
        do im = 1, basemodellist%Count()
          ! set model pointer
          num_model => GetNumericalModelFromList(basemodellist, im)
          if (num_model%name == export_model%modelname .and. &
              num_model%macronym == export_model%modeltype) then
            ! allocate and initialize nc export model
            call create_nc_export(export_model, num_model)
            exit
          end if
        end do
      end if
    end do
  end subroutine nc_export_create

end module NCExportCreateModule
