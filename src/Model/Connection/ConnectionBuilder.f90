module ConnectionBuilderModule
  use KindModule, only: I4B, LGP
  use SimModule, only: ustop
  use SimVariablesModule, only: iout
  use ListModule, only: ListType, isEqualIface
  use ListNodeModule, only: ListNodeType
  use BaseSolutionModule, only: BaseSolutionType
  use NumericalSolutionModule, only: NumericalSolutionType
  use BaseExchangeModule, only: BaseExchangeType, GetBaseExchangeFromList
  use DisConnExchangeModule, only: DisConnExchangeType, &
                                   GetDisConnExchangeFromList
  use NumericalModelModule, only: NumericalModelType
  use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                          cast_as_smc, &
                                          get_smc_from_list, &
                                          add_smc_to_list

  implicit none
  private

  type, public :: ConnectionBuilderType
    logical(LGP) :: dev_always_ifmod = .false. !< development option: force interface model on all exchanges
  contains
    procedure, pass(this) :: processSolution
    procedure, private, pass(this) :: processExchanges
    procedure, private, pass(this) :: setConnectionsToSolution
    procedure, private, pass(this) :: createModelConnectivity
  end type ConnectionBuilderType

contains

  !> @brief Process the exchanges in the solution into model connections
  !!
  !! This routine processes all exchanges in a solution and,
  !! when required, creates model connections of the proper
  !! type (GWF-GWF, GWT-GWT, ...) for a subset. It removes this
  !! subset of exchanges from the solution and replaces them with the
  !! created connections.
  !<
  subroutine processSolution(this, solution)
    class(ConnectionBuilderType) :: this !< the connection builder object
    class(BaseSolutionType), pointer :: solution !< the solution for which the exchanges are processed
    ! local
    class(NumericalSolutionType), pointer :: numSol
    type(ListType) :: newConnections

    ! we only deal with Num. Sol. here
    select type (solution)
    class is (NumericalSolutionType)
      numSol => solution
    class default
      return
    end select

    ! create the connections and add local exchanges
    call this%processExchanges(numSol%exchangelist, newConnections)
    if (newConnections%Count() == 0) then
      return
    end if

    write (iout, '(1x,a,i0,a,a)') 'Created ', newConnections%Count(), &
      ' model connections for solution ', trim(solution%name)

    ! create the topology of models participating in the interfaces
    call this%createModelConnectivity(newConnections)

    ! replace numerical exchanges in solution with connections
    call this%setConnectionsToSolution(newConnections, numSol)

    ! clean up local resources
    call newConnections%Clear(destroy=.false.)

  end subroutine processSolution

  !> @brief Create connections from exchanges
  !!
  !! If the configuration demands it, this will create connections,
  !! for the exchanges (one connection per exchange) add them to
  !! the global list, and return them as @param newConnections
  !<
  subroutine processExchanges(this, exchanges, newConnections)
    use ListsModule, only: baseconnectionlist, baseexchangelist
    class(ConnectionBuilderType) :: this !< the connection builder object
    type(ListType), pointer, intent(in) :: exchanges !< the list of exchanges to process
    type(ListType), intent(inout) :: newConnections !< the newly created connections
    ! local
    class(DisConnExchangeType), pointer :: conEx
    class(BaseExchangeType), pointer :: baseEx
    integer(I4B) :: iex, ibasex
    class(SpatialModelConnectionType), pointer :: modelConnection
    logical(LGP) :: isPeriodic

    do iex = 1, exchanges%Count()
      conEx => GetDisConnExchangeFromList(exchanges, iex)
      if (.not. associated(conEx)) then
        ! if it is not DisConnExchangeType, we can skip it
        cycle
      end if

      ! for now, if we have XT3D on the interface, we use a connection,
      ! (this will be more generic in the future)
      if (conEx%use_interface_model() .or. conEx%dev_ifmod_on &
          .or. this%dev_always_ifmod) then

        ! we should not get period connections here
        isPeriodic = (conEx%v_model1 == conEx%v_model2)
        if (isPeriodic) then
          write (*, *) 'Error (which should never happen): interface model '// &
            'does not support periodic boundary condition'
          call ustop()
        end if

        if (conEx%v_model1%is_local) then
          ! create model connection for model 1
          modelConnection => createModelConnection(conEx%model1, conEx)
          call add_smc_to_list(baseconnectionlist, modelConnection)
          call add_smc_to_list(newConnections, modelConnection)
        end if

        ! and for model 2
        if (conEx%v_model2%is_local) then
          modelConnection => createModelConnection(conEx%model2, conEx)
          call add_smc_to_list(baseconnectionlist, modelConnection)
          call add_smc_to_list(newConnections, modelConnection)
        end if

        ! remove this exchange from the base list, ownership
        ! now lies with the connection
        do ibasex = 1, baseexchangelist%Count()
          baseEx => GetBaseExchangeFromList(baseexchangelist, ibasex)
          if (conEx%id == baseEx%id) then
            call baseexchangelist%RemoveNode(ibasex, .false.)
            exit
          end if
        end do

      end if
    end do

  end subroutine processExchanges

  !> @brief Create a model connection of a given type
  !!
  !! This is a factory method to create the various types
  !! of model connections
  !<
  function createModelConnection(model, exchange) result(connection)
    use GwfGwfConnectionModule, only: GwfGwfConnectionType
    use GwtGwtConnectionModule, only: GwtGwtConnectionType
    use GweGweConnectionModule, only: GweGweConnectionType
    use GwfModule, only: GwfModelType

    class(NumericalModelType), pointer, intent(in) :: model !< the model for which the connection will be created
    class(DisConnExchangeType), pointer, intent(in) :: exchange !< the type of connection
    class(SpatialModelConnectionType), pointer :: connection !< the created connection

    ! different concrete connection types:
    class(GwfGwfConnectionType), pointer :: flowConnection => null()
    class(GwtGwtConnectionType), pointer :: transportConnection => null()
    class(GweGweConnectionType), pointer :: energyTransportConnection => null()

    connection => null()

    ! select on type of connection to create
    select case (exchange%typename)
    case ('GWF-GWF')
      allocate (GwfGwfConnectionType :: flowConnection)
      call flowConnection%construct(model, exchange)
      connection => flowConnection
      flowConnection => null()
    case ('GWT-GWT')
      allocate (GwtGwtConnectionType :: transportConnection)
      call transportConnection%construct(model, exchange)
      connection => transportConnection
      transportConnection => null()
    case ('GWE-GWE')
      allocate (GweGweConnectionType :: energyTransportConnection)
      call energyTransportConnection%construct(model, exchange)
      connection => energyTransportConnection
      energyTransportConnection => null()
    case default
      write (*, *) 'Error (which should never happen): '// &
        'undefined exchangetype found'
      call ustop()
    end select

  end function createModelConnection

  !> @brief Set connections to the solution
  !!
  !! This adds the connections to the solution and removes
  !! those exchanges which are replaced by a connection
  !<
  subroutine setConnectionsToSolution(this, connections, solution)
    class(ConnectionBuilderType) :: this !< the connection builder object
    type(ListType), intent(inout) :: connections !< the connections created for the solution
    class(NumericalSolutionType), pointer, intent(in) :: solution !< the solution to which the connections are set
    ! local
    type(ListType) :: keepList
    class(*), pointer :: exPtr, exPtr2, connPtr
    class(SpatialModelConnectionType), pointer :: conn
    integer(I4B) :: iex, iconn
    logical(LGP) :: keepExchange

    ! first add all exchanges not replaced by the connections to a list
    do iex = 1, solution%exchangelist%Count()
      exPtr => solution%exchangelist%GetItem(iex)
      ! will this exchange be replaced by a connection?
      keepExchange = .true.
      do iconn = 1, connections%Count()
        conn => get_smc_from_list(connections, iconn)
        exPtr2 => conn%prim_exchange
        if (associated(exPtr2, exPtr)) then
          ! if so, don't add it to the list
          keepExchange = .false.
          exit
        end if
      end do

      if (keepExchange) then
        call keepList%Add(exPtr)
      end if
    end do

    ! first add persisting exchanges
    call solution%exchangelist%Clear(destroy=.false.)
    do iex = 1, keepList%Count()
      exPtr => keepList%GetItem(iex)
      call solution%exchangelist%Add(exPtr)
    end do

    ! now add connections
    do iconn = 1, connections%Count()
      connPtr => connections%GetItem(iconn)
      call solution%exchangelist%Add(connPtr)
    end do

    ! clean up
    call keepList%Clear(destroy=.false.)

  end subroutine setConnectionsToSolution

  !> @brief Create connectivity of models which contribute to the interface
  !!
  !! This loops over all connections and creates a halo with all
  !! models from the numerical solution. The model halo will be used to
  !! extend the interface grid to include cells from models which are
  !< indirectly connected, through yet another exchange object.
  subroutine createModelConnectivity(this, connections)
    class(ConnectionBuilderType) :: this !< the connection builder object
    type(ListType), intent(inout) :: connections !< all connections that are created for this solution
    ! local
    integer(I4B) :: iconn
    class(SpatialModelConnectionType), pointer :: modelConn

    ! create halo for the model connections
    do iconn = 1, connections%Count()
      modelConn => get_smc_from_list(connections, iconn)
      call modelConn%createModelHalo()
    end do

  end subroutine createModelConnectivity

end module ConnectionBuilderModule
