# makefile created by pymake for the 'mf6' executable.


include ./makedefaults

# Define the source file directories
SOURCEDIR1=../src
SOURCEDIR2=../src/Distributed
SOURCEDIR3=../src/Exchange
SOURCEDIR4=../src/Idm
SOURCEDIR5=../src/Idm/selector
SOURCEDIR6=../src/Model
SOURCEDIR7=../src/Model/ChannelFlow
SOURCEDIR8=../src/Model/Connection
SOURCEDIR9=../src/Model/Discretization
SOURCEDIR10=../src/Model/Geometry
SOURCEDIR11=../src/Model/GroundWaterEnergy
SOURCEDIR12=../src/Model/GroundWaterFlow
SOURCEDIR13=../src/Model/GroundWaterFlow/submodules
SOURCEDIR14=../src/Model/GroundWaterTransport
SOURCEDIR15=../src/Model/ModelUtilities
SOURCEDIR16=../src/Model/OverlandFlow
SOURCEDIR17=../src/Model/ParticleTracking
SOURCEDIR18=../src/Model/SurfaceWaterFlow
SOURCEDIR19=../src/Model/TransportModel
SOURCEDIR20=../src/Solution
SOURCEDIR21=../src/Solution/LinearMethods
SOURCEDIR22=../src/Solution/PETSc
SOURCEDIR23=../src/Solution/ParticleTracker
SOURCEDIR24=../src/Timing
SOURCEDIR25=../src/Utilities
SOURCEDIR26=../src/Utilities/ArrayRead
SOURCEDIR27=../src/Utilities/Export
SOURCEDIR28=../src/Utilities/Idm
SOURCEDIR29=../src/Utilities/Idm/mf6blockfile
SOURCEDIR30=../src/Utilities/Idm/netcdf
SOURCEDIR31=../src/Utilities/Libraries
SOURCEDIR32=../src/Utilities/Libraries/blas
SOURCEDIR33=../src/Utilities/Libraries/daglib
SOURCEDIR34=../src/Utilities/Libraries/rcm
SOURCEDIR35=../src/Utilities/Libraries/sparsekit
SOURCEDIR36=../src/Utilities/Libraries/sparskit2
SOURCEDIR37=../src/Utilities/Matrix
SOURCEDIR38=../src/Utilities/Memory
SOURCEDIR39=../src/Utilities/Observation
SOURCEDIR40=../src/Utilities/OutputControl
SOURCEDIR41=../src/Utilities/Performance
SOURCEDIR42=../src/Utilities/TimeSeries
SOURCEDIR43=../src/Utilities/Vector

VPATH = \
${SOURCEDIR1} \
${SOURCEDIR2} \
${SOURCEDIR3} \
${SOURCEDIR4} \
${SOURCEDIR5} \
${SOURCEDIR6} \
${SOURCEDIR7} \
${SOURCEDIR8} \
${SOURCEDIR9} \
${SOURCEDIR10} \
${SOURCEDIR11} \
${SOURCEDIR12} \
${SOURCEDIR13} \
${SOURCEDIR14} \
${SOURCEDIR15} \
${SOURCEDIR16} \
${SOURCEDIR17} \
${SOURCEDIR18} \
${SOURCEDIR19} \
${SOURCEDIR20} \
${SOURCEDIR21} \
${SOURCEDIR22} \
${SOURCEDIR23} \
${SOURCEDIR24} \
${SOURCEDIR25} \
${SOURCEDIR26} \
${SOURCEDIR27} \
${SOURCEDIR28} \
${SOURCEDIR29} \
${SOURCEDIR30} \
${SOURCEDIR31} \
${SOURCEDIR32} \
${SOURCEDIR33} \
${SOURCEDIR34} \
${SOURCEDIR35} \
${SOURCEDIR36} \
${SOURCEDIR37} \
${SOURCEDIR38} \
${SOURCEDIR39} \
${SOURCEDIR40} \
${SOURCEDIR41} \
${SOURCEDIR42} \
${SOURCEDIR43} 

.SUFFIXES: .f90 .F90 .o

OBJECTS = \
$(OBJDIR)/kind.o \
$(OBJDIR)/Constants.o \
$(OBJDIR)/ErrorUtil.o \
$(OBJDIR)/CharString.o \
$(OBJDIR)/SimVariables.o \
$(OBJDIR)/ArrayHandlers.o \
$(OBJDIR)/Message.o \
$(OBJDIR)/defmacro.o \
$(OBJDIR)/compilerversion.o \
$(OBJDIR)/version.o \
$(OBJDIR)/Sim.o \
$(OBJDIR)/OpenSpec.o \
$(OBJDIR)/KeyValueNode.o \
$(OBJDIR)/Iterator.o \
$(OBJDIR)/InputOutput.o \
$(OBJDIR)/KeyValueListIterator.o \
$(OBJDIR)/TableTerm.o \
$(OBJDIR)/KeyValueList.o \
$(OBJDIR)/Table.o \
$(OBJDIR)/MemoryHelper.o \
$(OBJDIR)/ListNode.o \
$(OBJDIR)/PtrHashTableIterator.o \
$(OBJDIR)/Memory.o \
$(OBJDIR)/InputDefinition.o \
$(OBJDIR)/ListIterator.o \
$(OBJDIR)/PtrHashTable.o \
$(OBJDIR)/MemoryContainerIterator.o \
$(OBJDIR)/utl-ncfidm.o \
$(OBJDIR)/utl-hpcidm.o \
$(OBJDIR)/swf-zdgidm.o \
$(OBJDIR)/swf-stoidm.o \
$(OBJDIR)/swf-pcpidm.o \
$(OBJDIR)/swf-namidm.o \
$(OBJDIR)/swf-icidm.o \
$(OBJDIR)/swf-flwidm.o \
$(OBJDIR)/swf-evpidm.o \
$(OBJDIR)/swf-disv2didm.o \
$(OBJDIR)/swf-disv1didm.o \
$(OBJDIR)/swf-dis2didm.o \
$(OBJDIR)/swf-dfwidm.o \
$(OBJDIR)/swf-cxsidm.o \
$(OBJDIR)/swf-chdidm.o \
$(OBJDIR)/swf-cdbidm.o \
$(OBJDIR)/sim-tdisidm.o \
$(OBJDIR)/sim-namidm.o \
$(OBJDIR)/prt-namidm.o \
$(OBJDIR)/prt-mipidm.o \
$(OBJDIR)/prt-disvidm.o \
$(OBJDIR)/prt-disidm.o \
$(OBJDIR)/olf-zdgidm.o \
$(OBJDIR)/olf-stoidm.o \
$(OBJDIR)/olf-pcpidm.o \
$(OBJDIR)/olf-namidm.o \
$(OBJDIR)/olf-icidm.o \
$(OBJDIR)/olf-flwidm.o \
$(OBJDIR)/olf-evpidm.o \
$(OBJDIR)/olf-disv2didm.o \
$(OBJDIR)/olf-dis2didm.o \
$(OBJDIR)/olf-dfwidm.o \
$(OBJDIR)/olf-chdidm.o \
$(OBJDIR)/olf-cdbidm.o \
$(OBJDIR)/gwt-namidm.o \
$(OBJDIR)/gwt-mstidm.o \
$(OBJDIR)/gwt-istidm.o \
$(OBJDIR)/gwt-icidm.o \
$(OBJDIR)/gwt-dspidm.o \
$(OBJDIR)/gwt-disvidm.o \
$(OBJDIR)/gwt-disuidm.o \
$(OBJDIR)/gwt-disidm.o \
$(OBJDIR)/gwt-cncidm.o \
$(OBJDIR)/gwf-welidm.o \
$(OBJDIR)/gwf-stoidm.o \
$(OBJDIR)/gwf-rividm.o \
$(OBJDIR)/gwf-rchidm.o \
$(OBJDIR)/gwf-rchaidm.o \
$(OBJDIR)/gwf-npfidm.o \
$(OBJDIR)/gwf-namidm.o \
$(OBJDIR)/gwf-icidm.o \
$(OBJDIR)/gwf-ghbidm.o \
$(OBJDIR)/gwf-evtidm.o \
$(OBJDIR)/gwf-evtaidm.o \
$(OBJDIR)/gwf-drnidm.o \
$(OBJDIR)/gwf-disvidm.o \
$(OBJDIR)/gwf-disuidm.o \
$(OBJDIR)/gwf-disidm.o \
$(OBJDIR)/gwf-chdidm.o \
$(OBJDIR)/gwe-namidm.o \
$(OBJDIR)/gwe-icidm.o \
$(OBJDIR)/gwe-disvidm.o \
$(OBJDIR)/gwe-disuidm.o \
$(OBJDIR)/gwe-disidm.o \
$(OBJDIR)/gwe-ctpidm.o \
$(OBJDIR)/gwe-cndidm.o \
$(OBJDIR)/exg-olfgwfidm.o \
$(OBJDIR)/exg-gwtgwtidm.o \
$(OBJDIR)/exg-gwfprtidm.o \
$(OBJDIR)/exg-gwfgwtidm.o \
$(OBJDIR)/exg-gwfgwfidm.o \
$(OBJDIR)/exg-gwfgweidm.o \
$(OBJDIR)/exg-gwegweidm.o \
$(OBJDIR)/exg-chfgwfidm.o \
$(OBJDIR)/chf-zdgidm.o \
$(OBJDIR)/chf-stoidm.o \
$(OBJDIR)/chf-pcpidm.o \
$(OBJDIR)/chf-namidm.o \
$(OBJDIR)/chf-icidm.o \
$(OBJDIR)/chf-flwidm.o \
$(OBJDIR)/chf-evpidm.o \
$(OBJDIR)/chf-disv1didm.o \
$(OBJDIR)/chf-dfwidm.o \
$(OBJDIR)/chf-cxsidm.o \
$(OBJDIR)/chf-chdidm.o \
$(OBJDIR)/chf-cdbidm.o \
$(OBJDIR)/List.o \
$(OBJDIR)/LongLineReader.o \
$(OBJDIR)/DevFeature.o \
$(OBJDIR)/MemoryStore.o \
$(OBJDIR)/IdmUtlDfnSelector.o \
$(OBJDIR)/IdmSwfDfnSelector.o \
$(OBJDIR)/IdmSimDfnSelector.o \
$(OBJDIR)/IdmPrtDfnSelector.o \
$(OBJDIR)/IdmOlfDfnSelector.o \
$(OBJDIR)/IdmGwtDfnSelector.o \
$(OBJDIR)/IdmGwfDfnSelector.o \
$(OBJDIR)/IdmGweDfnSelector.o \
$(OBJDIR)/IdmExgDfnSelector.o \
$(OBJDIR)/IdmChfDfnSelector.o \
$(OBJDIR)/TimeSeriesRecord.o \
$(OBJDIR)/MathUtil.o \
$(OBJDIR)/BlockParser.o \
$(OBJDIR)/MemoryManager.o \
$(OBJDIR)/IdmDfnSelector.o \
$(OBJDIR)/TimeSeries.o \
$(OBJDIR)/SourceCommon.o \
$(OBJDIR)/MemoryManagerExt.o \
$(OBJDIR)/ats.o \
$(OBJDIR)/GeomUtil.o \
$(OBJDIR)/TimeSeriesLink.o \
$(OBJDIR)/TimeSeriesFileList.o \
$(OBJDIR)/tdis.o \
$(OBJDIR)/HashTable.o \
$(OBJDIR)/VectorBase.o \
$(OBJDIR)/Sparse.o \
$(OBJDIR)/DisvGeom.o \
$(OBJDIR)/ArrayReaders.o \
$(OBJDIR)/TimeSeriesManager.o \
$(OBJDIR)/SmoothingFunctions.o \
$(OBJDIR)/MatrixBase.o \
$(OBJDIR)/ListReader.o \
$(OBJDIR)/Connections.o \
$(OBJDIR)/TimeArray.o \
$(OBJDIR)/ObsOutput.o \
$(OBJDIR)/DiscretizationBase.o \
$(OBJDIR)/TimeArraySeries.o \
$(OBJDIR)/ObsOutputList.o \
$(OBJDIR)/Observe.o \
$(OBJDIR)/BudgetFileReader.o \
$(OBJDIR)/TimeArraySeriesLink.o \
$(OBJDIR)/ObsUtility.o \
$(OBJDIR)/ObsContainer.o \
$(OBJDIR)/Disv1dGeom.o \
$(OBJDIR)/BudgetTerm.o \
$(OBJDIR)/Budget.o \
$(OBJDIR)/TimeArraySeriesManager.o \
$(OBJDIR)/PackageMover.o \
$(OBJDIR)/Obs.o \
$(OBJDIR)/NumericalPackage.o \
$(OBJDIR)/Particle.o \
$(OBJDIR)/PackageBudget.o \
$(OBJDIR)/HeadFileReader.o \
$(OBJDIR)/GridFileReader.o \
$(OBJDIR)/Disv.o \
$(OBJDIR)/Disv2d.o \
$(OBJDIR)/Disv1d.o \
$(OBJDIR)/Disu.o \
$(OBJDIR)/Dis.o \
$(OBJDIR)/Dis2d.o \
$(OBJDIR)/BudgetObject.o \
$(OBJDIR)/BoundaryPackage.o \
$(OBJDIR)/CellDefn.o \
$(OBJDIR)/TrackFile.o \
$(OBJDIR)/sort.o \
$(OBJDIR)/FlowModelInterface.o \
$(OBJDIR)/Cell.o \
$(OBJDIR)/Subcell.o \
$(OBJDIR)/TrackControl.o \
$(OBJDIR)/TimeSelect.o \
$(OBJDIR)/prt-fmi.o \
$(OBJDIR)/TimeStepSelect.o \
$(OBJDIR)/SfrCrossSectionUtils.o \
$(OBJDIR)/TernarySolveTrack.o \
$(OBJDIR)/SubcellTri.o \
$(OBJDIR)/Method.o \
$(OBJDIR)/SubcellRect.o \
$(OBJDIR)/VirtualBase.o \
$(OBJDIR)/STLVecInt.o \
$(OBJDIR)/BaseModel.o \
$(OBJDIR)/PrintSaveManager.o \
$(OBJDIR)/SfrCrossSectionManager.o \
$(OBJDIR)/dag_module.o \
$(OBJDIR)/BoundaryPackageExt.o \
$(OBJDIR)/MethodSubcellTernary.o \
$(OBJDIR)/MethodSubcellPollock.o \
$(OBJDIR)/VirtualDataLists.o \
$(OBJDIR)/VirtualDataContainer.o \
$(OBJDIR)/SimStages.o \
$(OBJDIR)/NumericalModel.o \
$(OBJDIR)/OutputControlData.o \
$(OBJDIR)/gwf-ic.o \
$(OBJDIR)/UzfEtUtil.o \
$(OBJDIR)/Xt3dAlgorithm.o \
$(OBJDIR)/TvBase.o \
$(OBJDIR)/gwf-sfr.o \
$(OBJDIR)/gwf-riv.o \
$(OBJDIR)/gwf-maw.o \
$(OBJDIR)/mf6lists.o \
$(OBJDIR)/gwf-lak.o \
$(OBJDIR)/GwfVscInputData.o \
$(OBJDIR)/gwf-ghb.o \
$(OBJDIR)/gwf-drn.o \
$(OBJDIR)/IndexMap.o \
$(OBJDIR)/ArrayReaderBase.o \
$(OBJDIR)/MethodSubcellPool.o \
$(OBJDIR)/MethodCell.o \
$(OBJDIR)/CellPoly.o \
$(OBJDIR)/CellRectQuad.o \
$(OBJDIR)/CellRect.o \
$(OBJDIR)/SwfCxsUtils.o \
$(OBJDIR)/VirtualModel.o \
$(OBJDIR)/BaseExchange.o \
$(OBJDIR)/TspSpc.o \
$(OBJDIR)/tsp-fmi.o \
$(OBJDIR)/GweInputData.o \
$(OBJDIR)/OutputControl.o \
$(OBJDIR)/tsp-ic.o \
$(OBJDIR)/TspAdvOptions.o \
$(OBJDIR)/UzfCellGroup.o \
$(OBJDIR)/Xt3dInterface.o \
$(OBJDIR)/gwf-tvk.o \
$(OBJDIR)/SpdisWorkArray.o \
$(OBJDIR)/HGeoUtil.o \
$(OBJDIR)/gwf-vsc.o \
$(OBJDIR)/GwfNpfOptions.o \
$(OBJDIR)/GwfConductanceUtils.o \
$(OBJDIR)/InterfaceMap.o \
$(OBJDIR)/SeqVector.o \
$(OBJDIR)/ImsLinearSettings.o \
$(OBJDIR)/ConvergenceSummary.o \
$(OBJDIR)/Integer2dReader.o \
$(OBJDIR)/MethodCellTernary.o \
$(OBJDIR)/MethodCellPollockQuad.o \
$(OBJDIR)/MethodCellPollock.o \
$(OBJDIR)/MethodCellPassToBot.o \
$(OBJDIR)/VectorInterpolation.o \
$(OBJDIR)/swf-cxs.o \
$(OBJDIR)/CellWithNbrs.o \
$(OBJDIR)/NumericalExchange.o \
$(OBJDIR)/tsp-ssm.o \
$(OBJDIR)/tsp-oc.o \
$(OBJDIR)/tsp-obs.o \
$(OBJDIR)/tsp-mvt.o \
$(OBJDIR)/tsp-adv.o \
$(OBJDIR)/gwf-uzf.o \
$(OBJDIR)/tsp-apt.o \
$(OBJDIR)/gwt-mst.o \
$(OBJDIR)/GwtDspOptions.o \
$(OBJDIR)/gwf-tvs.o \
$(OBJDIR)/GwfStorageUtils.o \
$(OBJDIR)/Mover.o \
$(OBJDIR)/GwfMvrPeriodData.o \
$(OBJDIR)/ImsLinearMisc.o \
$(OBJDIR)/gwf-npf.o \
$(OBJDIR)/GwfBuyInputData.o \
$(OBJDIR)/GweCndOptions.o \
$(OBJDIR)/VirtualSolution.o \
$(OBJDIR)/SparseMatrix.o \
$(OBJDIR)/LinearSolverBase.o \
$(OBJDIR)/ImsReordering.o \
$(OBJDIR)/ModflowInput.o \
$(OBJDIR)/DefinitionSelect.o \
$(OBJDIR)/StructVector.o \
$(OBJDIR)/IdmLogger.o \
$(OBJDIR)/NCFileVars.o \
$(OBJDIR)/Integer1dReader.o \
$(OBJDIR)/Double2dReader.o \
$(OBJDIR)/Double1dReader.o \
$(OBJDIR)/MethodCellPool.o \
$(OBJDIR)/CellUtil.o \
$(OBJDIR)/swf-dfw.o \
$(OBJDIR)/swf-ic.o \
$(OBJDIR)/VirtualExchange.o \
$(OBJDIR)/GridSorting.o \
$(OBJDIR)/DisConnExchange.o \
$(OBJDIR)/CsrUtils.o \
$(OBJDIR)/tsp.o \
$(OBJDIR)/gwt-uzt.o \
$(OBJDIR)/gwt-src.o \
$(OBJDIR)/gwt-sft.o \
$(OBJDIR)/gwt-mwt.o \
$(OBJDIR)/gwt-lkt.o \
$(OBJDIR)/gwt-ist.o \
$(OBJDIR)/gwt-dsp.o \
$(OBJDIR)/gwt-cnc.o \
$(OBJDIR)/gwf-api.o \
$(OBJDIR)/gwf-wel.o \
$(OBJDIR)/gwf-rch.o \
$(OBJDIR)/gwf-sto.o \
$(OBJDIR)/gwf-oc.o \
$(OBJDIR)/gwf-obs.o \
$(OBJDIR)/gwf-mvr.o \
$(OBJDIR)/gwf-hfb.o \
$(OBJDIR)/gwf-csub.o \
$(OBJDIR)/gwf-buy.o \
$(OBJDIR)/GhostNode.o \
$(OBJDIR)/gwf-evt.o \
$(OBJDIR)/gwf-chd.o \
$(OBJDIR)/gwe-uze.o \
$(OBJDIR)/gwe-sfe.o \
$(OBJDIR)/gwe-mwe.o \
$(OBJDIR)/gwe-lke.o \
$(OBJDIR)/gwe-est.o \
$(OBJDIR)/gwe-esl.o \
$(OBJDIR)/gwe-ctp.o \
$(OBJDIR)/gwe-cnd.o \
$(OBJDIR)/RouterBase.o \
$(OBJDIR)/STLStackInt.o \
$(OBJDIR)/ImsLinearSolver.o \
$(OBJDIR)/ImsLinearBase.o \
$(OBJDIR)/DynamicPackageParams.o \
$(OBJDIR)/StructArray.o \
$(OBJDIR)/LoadNCInput.o \
$(OBJDIR)/LayeredArrayReader.o \
$(OBJDIR)/InputLoadType.o \
$(OBJDIR)/ReleaseSchedule.o \
$(OBJDIR)/MethodDisv.o \
$(OBJDIR)/MethodDis.o \
$(OBJDIR)/swf-zdg.o \
$(OBJDIR)/swf-sto.o \
$(OBJDIR)/swf-pcp.o \
$(OBJDIR)/swf-oc.o \
$(OBJDIR)/swf-obs.o \
$(OBJDIR)/swf-flw.o \
$(OBJDIR)/swf-evp.o \
$(OBJDIR)/swf-cdb.o \
$(OBJDIR)/GridConnection.o \
$(OBJDIR)/DistributedVariable.o \
$(OBJDIR)/gwt.o \
$(OBJDIR)/gwf.o \
$(OBJDIR)/VirtualGwfModel.o \
$(OBJDIR)/GwfExchangeMover.o \
$(OBJDIR)/gwe.o \
$(OBJDIR)/SerialRouter.o \
$(OBJDIR)/Timer.o \
$(OBJDIR)/Profiler.o \
$(OBJDIR)/LinearSolverFactory.o \
$(OBJDIR)/ImsLinear.o \
$(OBJDIR)/BaseSolution.o \
$(OBJDIR)/BoundInputContext.o \
$(OBJDIR)/LoadMf6File.o \
$(OBJDIR)/AsciiInputLoadType.o \
$(OBJDIR)/prt-prp.o \
$(OBJDIR)/prt-oc.o \
$(OBJDIR)/prt-mip.o \
$(OBJDIR)/MethodPool.o \
$(OBJDIR)/swf.o \
$(OBJDIR)/SpatialModelConnection.o \
$(OBJDIR)/GwtInterfaceModel.o \
$(OBJDIR)/exg-gwtgwt.o \
$(OBJDIR)/GwfInterfaceModel.o \
$(OBJDIR)/exg-gwfgwf.o \
$(OBJDIR)/GweInterfaceModel.o \
$(OBJDIR)/exg-gwegwe.o \
$(OBJDIR)/RouterFactory.o \
$(OBJDIR)/NumericalSolution.o \
$(OBJDIR)/MappedMemory.o \
$(OBJDIR)/NCModel.o \
$(OBJDIR)/Mf6FileStoInput.o \
$(OBJDIR)/Mf6FileListInput.o \
$(OBJDIR)/Mf6FileGridInput.o \
$(OBJDIR)/prt.o \
$(OBJDIR)/olf.o \
$(OBJDIR)/chf.o \
$(OBJDIR)/ExplicitSolution.o \
$(OBJDIR)/exg-swfgwf.o \
$(OBJDIR)/GwtGwtConnection.o \
$(OBJDIR)/GwfGwfConnection.o \
$(OBJDIR)/GweGweConnection.o \
$(OBJDIR)/VirtualDataManager.o \
$(OBJDIR)/Mapper.o \
$(OBJDIR)/ModelExport.o \
$(OBJDIR)/IdmMf6File.o \
$(OBJDIR)/ModelPackageInput.o \
$(OBJDIR)/VirtualGwtModel.o \
$(OBJDIR)/VirtualGwtExchange.o \
$(OBJDIR)/VirtualGwfExchange.o \
$(OBJDIR)/VirtualGweModel.o \
$(OBJDIR)/VirtualGweExchange.o \
$(OBJDIR)/SolutionGroup.o \
$(OBJDIR)/SolutionFactory.o \
$(OBJDIR)/exg-olfgwf.o \
$(OBJDIR)/exg-gwfprt.o \
$(OBJDIR)/exg-gwfgwt.o \
$(OBJDIR)/exg-gwfgwe.o \
$(OBJDIR)/DistributedSim.o \
$(OBJDIR)/exg-chfgwf.o \
$(OBJDIR)/RunControl.o \
$(OBJDIR)/SourceLoad.o \
$(OBJDIR)/ModelPackageInputs.o \
$(OBJDIR)/SimulationCreate.o \
$(OBJDIR)/RunControlFactory.o \
$(OBJDIR)/IdmLoad.o \
$(OBJDIR)/ConnectionBuilder.o \
$(OBJDIR)/comarg.o \
$(OBJDIR)/mf6core.o \
$(OBJDIR)/BaseGeometry.o \
$(OBJDIR)/mf6.o \
$(OBJDIR)/StringList.o \
$(OBJDIR)/MemorySetHandler.o \
$(OBJDIR)/ilut.o \
$(OBJDIR)/sparsekit.o \
$(OBJDIR)/rcm.o \
$(OBJDIR)/blas1_d.o \
$(OBJDIR)/Iunit.o \
$(OBJDIR)/GwfSfrCommon.o \
$(OBJDIR)/gwf-sfr-transient.o \
$(OBJDIR)/gwf-sfr-steady.o \
$(OBJDIR)/gwf-sfr-constant.o \
$(OBJDIR)/RectangularGeometry.o \
$(OBJDIR)/CircularGeometry.o \
$(OBJDIR)/ExplicitModel.o \
$(OBJDIR)/olf-disv1didm.o \
$(OBJDIR)/olf-cxsidm.o \
$(OBJDIR)/exg-swfgwfidm.o \
$(OBJDIR)/chf-disv2didm.o \
$(OBJDIR)/chf-dis2didm.o

# Define the objects that make up the program
$(PROGRAM) : $(OBJECTS)
	-$(FC) $(OPTLEVEL) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJDIR)/%.o : %.f90
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

$(OBJDIR)/%.o : %.F90
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

