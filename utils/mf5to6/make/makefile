# makefile created by pymake for the 'mf5to6' executable.


include ./makedefaults

# Define the source file directories
SOURCEDIR1=../src
SOURCEDIR2=../src/LGR
SOURCEDIR3=../src/MF2005
SOURCEDIR4=../src/NWT
SOURCEDIR5=../src/Preproc
SOURCEDIR6=../../../src/Utilities/Memory
SOURCEDIR7=../../../src/Utilities/TimeSeries
SOURCEDIR8=../../../src/Utilities

VPATH = \
${SOURCEDIR1} \
${SOURCEDIR2} \
${SOURCEDIR3} \
${SOURCEDIR4} \
${SOURCEDIR5} \
${SOURCEDIR6} \
${SOURCEDIR7} \
${SOURCEDIR8} 

.SUFFIXES: .f90 .F90 .f .o

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
$(OBJDIR)/PtrHashTableIterator.o \
$(OBJDIR)/Memory.o \
$(OBJDIR)/GlobalVariables.o \
$(OBJDIR)/Global.o \
$(OBJDIR)/PtrHashTable.o \
$(OBJDIR)/MemoryContainerIterator.o \
$(OBJDIR)/SimPHMF.o \
$(OBJDIR)/ListNode.o \
$(OBJDIR)/MemoryStore.o \
$(OBJDIR)/LongLineReader.o \
$(OBJDIR)/FeatureFlags.o \
$(OBJDIR)/Utilities.o \
$(OBJDIR)/ListIterator.o \
$(OBJDIR)/ConstantsPHMF.o \
$(OBJDIR)/MemoryManager.o \
$(OBJDIR)/GeomUtil.o \
$(OBJDIR)/BlockParser.o \
$(OBJDIR)/ArrayReadersMF5.o \
$(OBJDIR)/precutls.o \
$(OBJDIR)/MathUtil.o \
$(OBJDIR)/List.o \
$(OBJDIR)/GlobalVariablesPHMF.o \
$(OBJDIR)/DiscretizationBasePHMF.o \
$(OBJDIR)/File.o \
$(OBJDIR)/ObservePHMF.o \
$(OBJDIR)/Discretization3D.o \
$(OBJDIR)/CharacterContainer.o \
$(OBJDIR)/FileList.o \
$(OBJDIR)/ObsBlock.o \
$(OBJDIR)/LineList.o \
$(OBJDIR)/GlobalPHMF.o \
$(OBJDIR)/utl7.o \
$(OBJDIR)/StressPeriod.o \
$(OBJDIR)/Preproc.o \
$(OBJDIR)/obs2riv7.o \
$(OBJDIR)/obs2ghb7.o \
$(OBJDIR)/obs2drn7.o \
$(OBJDIR)/obs2chd7.o \
$(OBJDIR)/obs2bas7.o \
$(OBJDIR)/MultiLayerObsModule.o \
$(OBJDIR)/FileWriter.o \
$(OBJDIR)/ConverterCommon.o \
$(OBJDIR)/SimListVariables.o \
$(OBJDIR)/SfrDiversion.o \
$(OBJDIR)/Auxiliary.o \
$(OBJDIR)/TimeSeriesRecord.o \
$(OBJDIR)/ObsWriter.o \
$(OBJDIR)/Mover.o \
$(OBJDIR)/ModelPackage.o \
$(OBJDIR)/SfrReach.o \
$(OBJDIR)/GwfSfrModule.o \
$(OBJDIR)/ParamModule.o \
$(OBJDIR)/TimeSeries.o \
$(OBJDIR)/PackageWriter.o \
$(OBJDIR)/SfrSegment.o \
$(OBJDIR)/GwfSfrCheck.o \
$(OBJDIR)/GwfLpfModule.o \
$(OBJDIR)/GwfBcfModule.o \
$(OBJDIR)/GwfBasModule.o \
$(OBJDIR)/GwfLakModule.o \
$(OBJDIR)/GwfWelModule.o \
$(OBJDIR)/TdisWriter.o \
$(OBJDIR)/GwfChdModule.o \
$(OBJDIR)/ChdType.o \
$(OBJDIR)/MvrPackageWriter.o \
$(OBJDIR)/GwfSfrSubs.o \
$(OBJDIR)/GwfUzfModule_NWT.o \
$(OBJDIR)/GwfUpwModule.o \
$(OBJDIR)/GwfRivModule.o \
$(OBJDIR)/GwfRchModule.o \
$(OBJDIR)/GwfMnwModule.o \
$(OBJDIR)/LakeOutlet.o \
$(OBJDIR)/LakeConnection.o \
$(OBJDIR)/mach_mod.o \
$(OBJDIR)/gwf2hfb7_NWT.o \
$(OBJDIR)/GwfGhbModule.o \
$(OBJDIR)/gwf2wel7_NWT.o \
$(OBJDIR)/GwfFhbModule.o \
$(OBJDIR)/TdisVariables.o \
$(OBJDIR)/GwfChdSubs.o \
$(OBJDIR)/ChdObsWriter.o \
$(OBJDIR)/GwfEvtModule.o \
$(OBJDIR)/GwfEtsModule.o \
$(OBJDIR)/GwfDrnModule.o \
$(OBJDIR)/SfrPackageWriter.o \
$(OBJDIR)/GwfUzfSubs_NWT.o \
$(OBJDIR)/GwfRivSubs.o \
$(OBJDIR)/GwfRchSubs.o \
$(OBJDIR)/GwfBasOcSubs.o \
$(OBJDIR)/GwfMnwSubs.o \
$(OBJDIR)/LakeTributary.o \
$(OBJDIR)/Lake.o \
$(OBJDIR)/GwfLakSubs.o \
$(OBJDIR)/sip7.o \
$(OBJDIR)/pcgn2.o \
$(OBJDIR)/pcg7.o \
$(OBJDIR)/NWT1_module.o \
$(OBJDIR)/gmg7.o \
$(OBJDIR)/de47.o \
$(OBJDIR)/GwfHfbSubs_NWT.o \
$(OBJDIR)/GwfGhbSubs.o \
$(OBJDIR)/WelPackageWriter.o \
$(OBJDIR)/GwfFhbSubs.o \
$(OBJDIR)/ChdPackageWriter.o \
$(OBJDIR)/GwfEvtSubs.o \
$(OBJDIR)/GwfEtsSubs.o \
$(OBJDIR)/GwfDrnSubs.o \
$(OBJDIR)/UzfPackageWriter.o \
$(OBJDIR)/StoWriter.o \
$(OBJDIR)/RivPackageWriter.o \
$(OBJDIR)/RchPackageWriter.o \
$(OBJDIR)/OutputControlWriter.o \
$(OBJDIR)/NpfWriter.o \
$(OBJDIR)/MawPackageWriter.o \
$(OBJDIR)/GwfLgrModule.o \
$(OBJDIR)/LakPackageWriter.o \
$(OBJDIR)/ImsPackageWriter.o \
$(OBJDIR)/IcWriter.o \
$(OBJDIR)/HfbPackageWriter.o \
$(OBJDIR)/GwfResModule.o \
$(OBJDIR)/GhbPackageWriter.o \
$(OBJDIR)/FhbPackageWriter.o \
$(OBJDIR)/EvtPackageWriter.o \
$(OBJDIR)/DrnPackageWriter.o \
$(OBJDIR)/DisWriter.o \
$(OBJDIR)/Model.o \
$(OBJDIR)/GwfUpwSubs.o \
$(OBJDIR)/RivObsWriter.o \
$(OBJDIR)/GwfLpfSubs.o \
$(OBJDIR)/GwfLgrSubs.o \
$(OBJDIR)/GwfBcfSubs.o \
$(OBJDIR)/GwfBasSubs.o \
$(OBJDIR)/GhbObsWriter.o \
$(OBJDIR)/DrnObsWriter.o \
$(OBJDIR)/Connection.o \
$(OBJDIR)/ModelConverter.o \
$(OBJDIR)/Exchange.o \
$(OBJDIR)/ExchangeWriter.o \
$(OBJDIR)/SimFileWriter.o \
$(OBJDIR)/NWT1_xmdlib.o \
$(OBJDIR)/mf5to6.o \
$(OBJDIR)/SimVariablesPHMF.o \
$(OBJDIR)/NWT1_xmd.o \
$(OBJDIR)/NWT1_solver.o \
$(OBJDIR)/NWT1_ilupc_mod.o \
$(OBJDIR)/NWT1_gmres.o \
$(OBJDIR)/parutl7.o \
$(OBJDIR)/gwf2mnw2i7.o \
$(OBJDIR)/gwf2mnw17.o \
$(OBJDIR)/GwfResSubs.o

# Define the objects that make up the program
$(PROGRAM) : $(OBJECTS)
	-$(FC) $(OPTLEVEL) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJDIR)/%.o : %.f90
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

$(OBJDIR)/%.o : %.F90
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

$(OBJDIR)/%.o : %.f
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

