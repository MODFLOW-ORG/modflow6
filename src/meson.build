external_libraries = files(
    'Utilities' / 'Libraries' / 'daglib' / 'dag_module.f90',
    'Utilities' / 'Libraries' / 'rcm' / 'rcm.f90',
    'Utilities' / 'Libraries' / 'sparsekit' / 'sparsekit.f90',
    'Utilities' / 'Libraries' / 'blas' / 'blas1_d.f90',
    'Utilities' / 'Libraries' / 'sparskit2' / 'ilut.f90',
)

modflow_sources = files(
    'Distributed' / 'DistributedSim.f90',
    'Distributed' / 'IndexMap.f90',
    'Distributed' / 'InterfaceMap.f90',
    'Distributed' / 'Mapper.f90',
    'Distributed' / 'MappedMemory.f90',
    'Distributed' / 'RouterBase.f90',
    'Distributed' / 'RouterFactory.F90',
    'Distributed' / 'SerialRouter.f90',
    'Distributed' / 'VirtualBase.f90',
    'Distributed' / 'VirtualDataContainer.f90',
    'Distributed' / 'VirtualDataLists.f90',
    'Distributed' / 'VirtualDataManager.f90',
    'Distributed' / 'VirtualExchange.f90',    
    'Distributed' / 'VirtualGweExchange.f90',
    'Distributed' / 'VirtualGwfExchange.f90',
    'Distributed' / 'VirtualGwtExchange.f90',
    'Distributed' / 'VirtualModel.f90',
    'Distributed' / 'VirtualGweModel.f90',
    'Distributed' / 'VirtualGwfModel.f90',
    'Distributed' / 'VirtualGwtModel.f90',
    'Distributed' / 'VirtualSolution.f90',
    'Exchange' / 'BaseExchange.f90',
    'Exchange' / 'DisConnExchange.f90',
    'Exchange' / 'GhostNode.f90',
    'Exchange' / 'GwfExchangeMover.f90',
    'Exchange' / 'NumericalExchange.f90',
    'Exchange' / 'exg-chfgwf.f90',
    'Exchange' / 'exg-gwfgwf.f90',
    'Exchange' / 'exg-gwfgwt.f90',
    'Exchange' / 'exg-gwtgwt.f90',
    'Exchange' / 'exg-gwfgwe.f90',
    'Exchange' / 'exg-gwegwe.f90',
    'Exchange' / 'exg-gwfprt.f90',
    'Exchange' / 'exg-olfgwf.f90',
    'Exchange' / 'exg-swfgwf.f90',
    'Idm' / 'exg-chfgwfidm.f90',
    'Idm' / 'exg-gwfgwfidm.f90',
    'Idm' / 'exg-gwfgwtidm.f90',
    'Idm' / 'exg-gwtgwtidm.f90',
    'Idm' / 'exg-gwfgweidm.f90',
    'Idm' / 'exg-gwegweidm.f90',
    'Idm' / 'exg-gwfprtidm.f90',
    'Idm' / 'exg-olfgwfidm.f90',
    'Idm' / 'gwe-cndidm.f90',
    'Idm' / 'gwe-ctpidm.f90',
    'Idm' / 'gwe-disidm.f90',
    'Idm' / 'gwe-disuidm.f90',
    'Idm' / 'gwe-disvidm.f90',
    'Idm' / 'gwe-icidm.f90',
    'Idm' / 'gwe-namidm.f90',
    'Idm' / 'gwf-chdidm.f90',
    'Idm' / 'gwf-disidm.f90',
    'Idm' / 'gwf-disuidm.f90',
    'Idm' / 'gwf-disvidm.f90',
    'Idm' / 'gwf-drnidm.f90',
    'Idm' / 'gwf-evtidm.f90',
    'Idm' / 'gwf-evtaidm.f90',
    'Idm' / 'gwf-ghbidm.f90',
    'Idm' / 'gwf-icidm.f90',
    'Idm' / 'gwf-namidm.f90',
    'Idm' / 'gwf-npfidm.f90',
    'Idm' / 'gwf-rchidm.f90',
    'Idm' / 'gwf-rchaidm.f90',
    'Idm' / 'gwf-rividm.f90',
    'Idm' / 'gwf-stoidm.f90',
    'Idm' / 'gwf-welidm.f90',
    'Idm' / 'gwt-cncidm.f90',
    'Idm' / 'gwt-disidm.f90',
    'Idm' / 'gwt-disuidm.f90',
    'Idm' / 'gwt-disvidm.f90',
    'Idm' / 'gwt-dspidm.f90',
    'Idm' / 'gwt-icidm.f90',
    'Idm' / 'gwt-namidm.f90',
    'Idm' / 'swf-disv1didm.f90',
    'Idm' / 'swf-dis2didm.f90',
    'Idm' / 'swf-disv2didm.f90',
    'Idm' / 'swf-namidm.f90',
    'Idm' / 'swf-cxsidm.f90',
    'Idm' / 'swf-dfwidm.f90',
    'Idm' / 'swf-icidm.f90',
    'Idm' / 'swf-flwidm.f90',
    'Idm' / 'swf-chdidm.f90',
    'Idm' / 'swf-cdbidm.f90',
    'Idm' / 'swf-evpidm.f90',
    'Idm' / 'swf-pcpidm.f90',
    'Idm' / 'swf-stoidm.f90',
    'Idm' / 'swf-zdgidm.f90',
    'Idm' / 'chf-disv1didm.f90',
    'Idm' / 'chf-dis2didm.f90',
    'Idm' / 'chf-disv2didm.f90',
    'Idm' / 'chf-namidm.f90',
    'Idm' / 'chf-cxsidm.f90',
    'Idm' / 'chf-dfwidm.f90',
    'Idm' / 'chf-icidm.f90',
    'Idm' / 'chf-flwidm.f90',
    'Idm' / 'chf-chdidm.f90',
    'Idm' / 'chf-cdbidm.f90',
    'Idm' / 'chf-evpidm.f90',
    'Idm' / 'chf-pcpidm.f90',
    'Idm' / 'chf-stoidm.f90',
    'Idm' / 'chf-zdgidm.f90',
    'Idm' / 'olf-disv1didm.f90',
    'Idm' / 'olf-dis2didm.f90',
    'Idm' / 'olf-disv2didm.f90',
    'Idm' / 'olf-namidm.f90',
    'Idm' / 'olf-cxsidm.f90',
    'Idm' / 'olf-dfwidm.f90',
    'Idm' / 'olf-icidm.f90',
    'Idm' / 'olf-flwidm.f90',
    'Idm' / 'olf-chdidm.f90',
    'Idm' / 'olf-cdbidm.f90',
    'Idm' / 'olf-evpidm.f90',
    'Idm' / 'olf-pcpidm.f90',
    'Idm' / 'olf-stoidm.f90',
    'Idm' / 'olf-zdgidm.f90',
    'Idm' / 'sim-namidm.f90',
    'Idm' / 'sim-tdisidm.f90',
    'Idm' / 'prt-namidm.f90',
    'Idm' / 'prt-disidm.f90',
    'Idm' / 'prt-disvidm.f90',
    'Idm' / 'prt-mipidm.f90',
    'Idm' / 'utl-hpcidm.f90',
    'Idm' / 'utl-ncfidm.f90',
    'Idm' / 'selector' / 'IdmDfnSelector.f90',
    'Idm' / 'selector' / 'IdmExgDfnSelector.f90',
    'Idm' / 'selector' / 'IdmGweDfnSelector.f90',
    'Idm' / 'selector' / 'IdmGwfDfnSelector.f90',
    'Idm' / 'selector' / 'IdmGwtDfnSelector.f90',
    'Idm' / 'selector' / 'IdmChfDfnSelector.f90',
    'Idm' / 'selector' / 'IdmSwfDfnSelector.f90',
    'Idm' / 'selector' / 'IdmOlfDfnSelector.f90',
    'Idm' / 'selector' / 'IdmPrtDfnSelector.f90',
    'Idm' / 'selector' / 'IdmSimDfnSelector.f90',
    'Idm' / 'selector' / 'IdmUtlDfnSelector.f90',
    'Model' / 'ChannelFlow' / 'chf.f90',
    'Model' / 'Connection' / 'ConnectionBuilder.f90',
    'Model' / 'Connection' / 'CellWithNbrs.f90',
    'Model' / 'Connection' / 'CsrUtils.f90',
    'Model' / 'Connection' / 'GridConnection.f90',
    'Model' / 'Connection' / 'GridSorting.f90',
    'Model' / 'Connection' / 'GweGweConnection.f90',
    'Model' / 'Connection' / 'GweInterfaceModel.f90',
    'Model' / 'Connection' / 'GwfGwfConnection.f90',
    'Model' / 'Connection' / 'GwtGwtConnection.f90',
    'Model' / 'Connection' / 'GwfInterfaceModel.f90',
    'Model' / 'Connection' / 'GwtInterfaceModel.f90',
    'Model' / 'Connection' / 'SpatialModelConnection.f90',
    'Model' / 'Connection' / 'DistributedVariable.f90',
    'Model' / 'Discretization' / 'Dis.f90',
    'Model' / 'Discretization' / 'Dis2d.f90',
    'Model' / 'Discretization' / 'Disu.f90',
    'Model' / 'Discretization' / 'Disv.f90',
    'Model' / 'Discretization' / 'Disv1d.f90',
    'Model' / 'Discretization' / 'Disv2d.f90',
    'Model' / 'Geometry' / 'BaseGeometry.f90',
    'Model' / 'Geometry' / 'CircularGeometry.f90',
    'Model' / 'Geometry' / 'RectangularGeometry.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-cnd.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-ctp.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-esl.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-est.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-lke.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-mwe.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-sfe.f90',
    'Model' / 'GroundWaterEnergy' / 'gwe-uze.f90',
    'Model' / 'GroundWaterFlow' / 'gwf.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-api.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-buy.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-chd.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-csub.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-drn.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-evt.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-ghb.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-hfb.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-ic.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-lak.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-maw.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-mvr.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-npf.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-obs.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-oc.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-rch.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-riv.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-sfr.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-sto.f90',
    'Model' / 'GroundWaterFlow' / 'TvBase.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-tvk.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-tvs.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-uzf.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-vsc.f90',
    'Model' / 'GroundWaterFlow' / 'gwf-wel.f90',
    'Model' / 'GroundWaterFlow' / 'submodules' / 'gwf-sfr-steady.f90',
    'Model' / 'GroundWaterFlow' / 'submodules' / 'gwf-sfr-transient.f90',
    'Model' / 'GroundWaterFlow' / 'submodules' / 'gwf-sfr-constant.f90',
    'Model' / 'GroundWaterTransport' / 'gwt.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-cnc.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-dsp.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-ist.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-lkt.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-mst.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-mwt.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-sft.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-src.f90',
    'Model' / 'GroundWaterTransport' / 'gwt-uzt.f90',
    'Model' / 'OverlandFlow' / 'olf.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-cxs.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-dfw.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-ic.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-obs.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-oc.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-flw.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-evp.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-pcp.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-sto.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-cdb.f90',
    'Model' / 'SurfaceWaterFlow' / 'swf-zdg.f90',
    'Model' / 'ParticleTracking' / 'prt.f90',
    'Model' / 'ParticleTracking' / 'prt-fmi.f90',
    'Model' / 'ParticleTracking' / 'prt-mip.f90',
    'Model' / 'ParticleTracking' / 'prt-oc.f90',
    'Model' / 'ParticleTracking' / 'prt-prp.f90',
    'Model' / 'ModelUtilities' / 'BoundaryPackage.f90',
    'Model' / 'ModelUtilities' / 'BoundaryPackageExt.f90',
    'Model' / 'ModelUtilities' / 'GweCndOptions.f90',
    'Model' / 'ModelUtilities' / 'Connections.f90',
    'Model' / 'ModelUtilities' / 'DiscretizationBase.f90',
    'Model' / 'ModelUtilities' / 'DisvGeom.f90',
    'Model' / 'ModelUtilities' / 'Disv1dGeom.f90',
    'Model' / 'ModelUtilities' / 'FlowModelInterface.f90',
    'Model' / 'ModelUtilities' / 'GweInputData.f90',
    'Model' / 'ModelUtilities' / 'GwfBuyInputData.f90',
    'Model' / 'ModelUtilities' / 'GwfConductanceUtils.f90',
    'Model' / 'ModelUtilities' / 'GwfMvrPeriodData.f90',
    'Model' / 'ModelUtilities' / 'GwfNpfOptions.f90',
    'Model' / 'ModelUtilities' / 'GwfSfrCommon.f90',
    'Model' / 'ModelUtilities' / 'GwfStorageUtils.f90',
    'Model' / 'ModelUtilities' / 'GwfVscInputData.f90',
    'Model' / 'ModelUtilities' / 'GwtDspOptions.f90',
    'Model' / 'ModelUtilities' / 'ModelPackageInput.f90',
    'Model' / 'ModelUtilities' / 'Mover.f90',
    'Model' / 'ModelUtilities' / 'PackageMover.f90',
    'Model' / 'ModelUtilities' / 'ReleaseSchedule.f90',
    'Model' / 'ModelUtilities' / 'SfrCrossSectionManager.f90',
    'Model' / 'ModelUtilities' / 'SfrCrossSectionUtils.f90',
    'Model' / 'ModelUtilities' / 'SpdisWorkArray.f90',
    'Model' / 'ModelUtilities' / 'SwfCxsUtils.f90',
    'Model' / 'ModelUtilities' / 'TspAdvOptions.f90',
    'Model' / 'ModelUtilities' / 'TspSpc.f90',
    'Model' / 'ModelUtilities' / 'TimeSelect.f90',
    'Model' / 'ModelUtilities' / 'TimeStepSelect.f90',
    'Model' / 'ModelUtilities' / 'TrackFile.f90',
    'Model' / 'ModelUtilities' / 'TrackControl.f90',
    'Model' / 'ModelUtilities' / 'UzfCellGroup.f90',
    'Model' / 'ModelUtilities' / 'UzfEtUtil.f90',
    'Model' / 'ModelUtilities' / 'VectorInterpolation.f90',
    'Model' / 'ModelUtilities' / 'Xt3dAlgorithm.f90',
    'Model' / 'ModelUtilities' / 'Xt3dInterface.f90',
    'Model' / 'TransportModel' / 'tsp.f90',
    'Model' / 'TransportModel' / 'tsp-adv.f90',
    'Model' / 'TransportModel' / 'tsp-apt.f90',
    'Model' / 'TransportModel' / 'tsp-fmi.f90',
    'Model' / 'TransportModel' / 'tsp-ic.f90',
    'Model' / 'TransportModel' / 'tsp-obs.f90',
    'Model' / 'TransportModel' / 'tsp-oc.f90',
    'Model' / 'TransportModel' / 'tsp-mvt.f90',
    'Model' / 'TransportModel' / 'tsp-ssm.f90',
    'Model' / 'BaseModel.f90',
    'Model' / 'ExplicitModel.f90',
    'Model' / 'NumericalModel.f90',
    'Model' / 'NumericalPackage.f90',
    'Solution' / 'ConvergenceSummary.f90',
    'Solution' / 'LinearMethods' / 'ImsLinearSettings.f90',
    'Solution' / 'LinearMethods' / 'ImsLinearBase.f90',
    'Solution' / 'LinearMethods' / 'ImsLinear.f90',
    'Solution' / 'LinearMethods' / 'ImsReordering.f90',
    'Solution' / 'LinearMethods' / 'ImsLinearMisc.f90',
    'Solution' / 'LinearMethods' / 'ImsLinearSolver.f90',
    'Solution' / 'BaseSolution.f90',
    'Solution' / 'ExplicitSolution.f90',
    'Solution' / 'LinearSolverBase.f90',
    'Solution' / 'LinearSolverFactory.F90',
    'Solution' / 'NumericalSolution.f90',
    'Solution' / 'SolutionFactory.F90',
    'Solution' / 'SolutionGroup.f90',
    'Solution' / 'ParticleTracker' / 'Cell.f90',
    'Solution' / 'ParticleTracker' / 'CellDefn.f90',
    'Solution' / 'ParticleTracker' / 'CellPoly.f90',
    'Solution' / 'ParticleTracker' / 'CellRect.f90',
    'Solution' / 'ParticleTracker' / 'CellRectQuad.f90',
    'Solution' / 'ParticleTracker' / 'CellUtil.f90',
    'Solution' / 'ParticleTracker' / 'Method.f90',
    'Solution' / 'ParticleTracker' / 'MethodCell.f90',
    'Solution' / 'ParticleTracker' / 'MethodCellPollock.f90',
    'Solution' / 'ParticleTracker' / 'MethodCellPollockQuad.f90',
    'Solution' / 'ParticleTracker' / 'MethodCellPool.f90',
    'Solution' / 'ParticleTracker' / 'MethodCellTernary.f90',
    'Solution' / 'ParticleTracker' / 'MethodDis.f90',
    'Solution' / 'ParticleTracker' / 'MethodDisv.f90',
    'Solution' / 'ParticleTracker' / 'MethodPool.f90',
    'Solution' / 'ParticleTracker' / 'MethodCellPassToBot.f90',
    'Solution' / 'ParticleTracker' / 'MethodSubcellPollock.f90',
    'Solution' / 'ParticleTracker' / 'MethodSubcellPool.f90',
    'Solution' / 'ParticleTracker' / 'MethodSubcellTernary.f90',
    'Solution' / 'ParticleTracker' / 'Particle.f90',
    'Solution' / 'ParticleTracker' / 'Subcell.f90',
    'Solution' / 'ParticleTracker' / 'SubcellRect.f90',
    'Solution' / 'ParticleTracker' / 'SubcellTri.f90',
    'Solution' / 'ParticleTracker' / 'TernarySolveTrack.f90',
    'Timing' / 'ats.f90',
    'Timing' / 'tdis.f90',
    'Utilities' / 'ArrayRead' / 'ArrayReaderBase.f90',
    'Utilities' / 'ArrayRead' / 'Double1dReader.f90',
    'Utilities' / 'ArrayRead' / 'Double2dReader.f90',
    'Utilities' / 'ArrayRead' / 'Integer1dReader.f90',
    'Utilities' / 'ArrayRead' / 'Integer2dReader.f90',
    'Utilities' / 'ArrayRead' / 'LayeredArrayReader.f90',
    'Utilities' / 'Export' / 'ModelExport.f90',
    'Utilities' / 'Export' / 'NCModel.f90',
    'Utilities' / 'Idm' / 'BoundInputContext.f90',
    'Utilities' / 'Idm' / 'DefinitionSelect.f90',
    'Utilities' / 'Idm' / 'DynamicPackageParams.f90',
    'Utilities' / 'Idm' / 'IdmLoad.f90',
    'Utilities' / 'Idm' / 'IdmLogger.f90',
    'Utilities' / 'Idm' / 'InputDefinition.f90',
    'Utilities' / 'Idm' / 'InputLoadType.f90',
    'Utilities' / 'Idm' / 'ModelPackageInputs.f90',
    'Utilities' / 'Idm' / 'ModflowInput.f90',
    'Utilities' / 'Idm' / 'SourceCommon.f90',
    'Utilities' / 'Idm' / 'SourceLoad.F90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'AsciiInputLoadType.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'IdmMf6File.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'LoadMf6File.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'Mf6FileGridInput.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'Mf6FileListInput.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'Mf6FileStoInput.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'LoadNCInput.F90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'StructArray.f90',
    'Utilities' / 'Idm' / 'mf6blockfile' / 'StructVector.f90',
    'Utilities' / 'Idm' / 'netcdf' / 'NCFileVars.f90',
    'Utilities' / 'Matrix' / 'MatrixBase.f90',
    'Utilities' / 'Matrix' / 'SparseMatrix.f90',
    'Utilities' / 'Memory' / 'Memory.F90',
    'Utilities' / 'Memory' / 'MemoryHelper.f90',
    'Utilities' / 'Memory' / 'MemoryContainerIterator.f90',
    'Utilities' / 'Memory' / 'MemoryStore.f90',
    'Utilities' / 'Memory' / 'MemoryManager.f90',
    'Utilities' / 'Memory' / 'MemoryManagerExt.f90',
    'Utilities' / 'Memory' / 'MemorySetHandler.f90',
    'Utilities' / 'Observation' / 'Obs.f90',
    'Utilities' / 'Observation' / 'ObsContainer.f90',
    'Utilities' / 'Observation' / 'Observe.f90',
    'Utilities' / 'Observation' / 'ObsOutput.f90',
    'Utilities' / 'Observation' / 'ObsOutputList.f90',
    'Utilities' / 'Observation' / 'ObsUtility.f90',
    'Utilities' / 'OutputControl' / 'OutputControl.f90',
    'Utilities' / 'OutputControl' / 'OutputControlData.f90',
    'Utilities' / 'OutputControl' / 'PrintSaveManager.f90',
    'Utilities' / 'Performance' / 'Profiler.f90',
    'Utilities' / 'TimeSeries' / 'TimeArray.f90',
    'Utilities' / 'TimeSeries' / 'TimeArraySeries.f90',
    'Utilities' / 'TimeSeries' / 'TimeArraySeriesLink.f90',
    'Utilities' / 'TimeSeries' / 'TimeArraySeriesManager.f90',
    'Utilities' / 'TimeSeries' / 'TimeSeries.f90',
    'Utilities' / 'TimeSeries' / 'TimeSeriesFileList.f90',
    'Utilities' / 'TimeSeries' / 'TimeSeriesLink.f90',
    'Utilities' / 'TimeSeries' / 'TimeSeriesManager.f90',
    'Utilities' / 'TimeSeries' / 'TimeSeriesRecord.f90',
    'Utilities' / 'Vector' / 'VectorBase.f90',
    'Utilities' / 'Vector' / 'SeqVector.f90',
    'Utilities' / 'ArrayHandlers.f90',
    'Utilities' / 'ArrayReaders.f90',
    'Utilities' / 'BlockParser.f90',
    'Utilities' / 'Budget.f90',
    'Utilities' / 'BudgetFileReader.f90',
    'Utilities' / 'BudgetObject.f90',
    'Utilities' / 'BudgetTerm.f90',
    'Utilities' / 'CharString.f90',
    'Utilities' / 'comarg.f90',
    'Utilities' / 'compilerversion.F90',
    'Utilities' / 'CharString.f90',
    'Utilities' / 'Constants.f90',
    'Utilities' / 'defmacro.F90',
    'Utilities' / 'DevFeature.f90',
    'Utilities' / 'ErrorUtil.f90',
    'Utilities' / 'GeomUtil.f90',
    'Utilities' / 'GridFileReader.f90',
    'Utilities' / 'HashTable.f90',
    'Utilities' / 'HeadFileReader.f90',
    'Utilities' / 'HGeoUtil.f90',
    'Utilities' / 'Iterator.f90',
    'Utilities' / 'InputOutput.f90',
    'Utilities' / 'Iunit.f90',
    'Utilities' / 'kind.f90',
    'Utilities' / 'KeyValueList.f90',
    'Utilities' / 'KeyValueListIterator.f90',
    'Utilities' / 'KeyValueNode.f90',
    'Utilities' / 'List.f90',
    'Utilities' / 'ListIterator.f90',
    'Utilities' / 'ListNode.f90',
    'Utilities' / 'ListReader.f90',
    'Utilities' / 'LongLineReader.f90',
    'Utilities' / 'MathUtil.f90',
    'Utilities' / 'Message.f90',
    'Utilities' / 'OpenSpec.f90',
    'Utilities' / 'PackageBudget.f90',
    'Utilities' / 'PtrHashTable.f90',
    'Utilities' / 'PtrHashTableIterator.f90',
    'Utilities' / 'Sim.f90',
    'Utilities' / 'SimStages.f90',
    'Utilities' / 'SimVariables.f90',
    'Utilities' / 'SmoothingFunctions.f90',
    'Utilities' / 'sort.f90',
    'Utilities' / 'Sparse.f90',
    'Utilities' / 'STLStackInt.f90',
    'Utilities' / 'STLVecInt.f90',
    'Utilities' / 'StringList.f90',
    'Utilities' / 'Table.f90',
    'Utilities' / 'TableTerm.f90',
    'Utilities' / 'Timer.f90',
    'Utilities' / 'version.f90',
    'mf6core.f90',
    'mf6lists.f90',
    'SimulationCreate.f90',
    'RunControl.f90',
    'RunControlFactory.F90'
)

modflow_petsc_sources = files(
  'Utilities' / 'Vector' / 'PetscVector.F90',
  'Utilities' / 'Matrix' / 'PetscMatrix.F90',
  'Solution' / 'PETSc' / 'PetscSolver.F90',
  'Solution' / 'PETSc' / 'PetscConvergence.F90',
  'Solution' / 'PETSc' / 'PetscImsPreconditioner.F90'
)

modflow_mpi_sources = files(
  'Distributed' / 'MpiMessageBuilder.f90',
  'Distributed' / 'MpiMessageCache.f90',
  'Distributed' / 'MpiUnitCache.f90',
  'Distributed' / 'MpiRouter.f90',
  'Distributed' / 'MpiRunControl.F90',
  'Distributed' / 'MpiWorld.f90',
  'Solution' / 'ParallelSolution.f90'
)

modflow_netcdf_sources = files(
  'Utilities' / 'Export' / 'DisNCMesh.f90',
  'Utilities' / 'Export' / 'DisNCStructured.f90',
  'Utilities' / 'Export' / 'DisvNCMesh.f90',
  'Utilities' / 'Export' / 'MeshNCModel.f90',
  'Utilities' / 'Export' / 'NCExportCreate.f90',
  'Utilities' / 'Idm' / 'netcdf' / 'NCArrayReader.f90',
  'Utilities' / 'Idm' / 'netcdf' / 'NCContextBuild.f90',
  'Utilities' / 'Idm' / 'netcdf' / 'NetCDFCommon.f90'
)

if with_petsc
  modflow_sources += modflow_petsc_sources  
endif

if with_mpi
  modflow_sources += modflow_mpi_sources
endif

if with_netcdf
  modflow_sources += modflow_netcdf_sources
endif

mf6_external = static_library('mf6_external', external_libraries)

message('MODFLOW 6 executable name: ' + buildname)

if build_machine.system() == 'windows' and (with_petsc or with_netcdf)
  if is_extended_build
    incdir = include_directories(extended_incdir)
  elif with_petsc
    incdir = petsc_incdir
  elif with_netcdf
    incdir = netcdf_incdir
  endif
  mf6core = static_library('mf6core',
                           modflow_sources,
                           dependencies:  dependencies,
                           link_with: [mf6_external],
                           include_directories: incdir)
  mf6exe = executable(buildname,
                      'mf6.f90',
                      link_with: [mf6core],
                      dependencies: dependencies,
                      install: true)
else
  mf6core = static_library('mf6core', 
                           modflow_sources,  
                           dependencies: dependencies,
                           link_with: [mf6_external])
  mf6exe = executable(buildname, 
                      'mf6.f90', 
                      link_with: [mf6core],  
                      dependencies: dependencies, 
                      install: true)
endif
