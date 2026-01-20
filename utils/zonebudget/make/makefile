# makefile created by pymake for the 'zbud6' executable.


include ./makedefaults

# Define the source file directories
SOURCEDIR1=../src
SOURCEDIR2=../../../src/Utilities

VPATH = \
${SOURCEDIR1} \
${SOURCEDIR2} 

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
$(OBJDIR)/InputOutput.o \
$(OBJDIR)/LongLineReader.o \
$(OBJDIR)/FeatureFlags.o \
$(OBJDIR)/sort.o \
$(OBJDIR)/BlockParser.o \
$(OBJDIR)/ArrayReaders.o \
$(OBJDIR)/zone.o \
$(OBJDIR)/Budget.o \
$(OBJDIR)/HashTable.o \
$(OBJDIR)/zoneoutput.o \
$(OBJDIR)/GridFileReader.o \
$(OBJDIR)/budgetdata.o \
$(OBJDIR)/MathUtil.o \
$(OBJDIR)/GeomUtil.o \
$(OBJDIR)/zbud6.o

# Define the objects that make up the program
$(PROGRAM) : $(OBJECTS)
	-$(FC) $(OPTLEVEL) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJDIR)/%.o : %.f90
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

$(OBJDIR)/%.o : %.F90
	@mkdir -p $(@D)
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ $(INCSWITCH) $(MODSWITCH)

