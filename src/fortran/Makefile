CODEVLIBDIR = C:/CODEV114/umr/lib

BUILDDIR = fortran/build

BINDIR = $(BUILDDIR)/bin
INCLUDEDIR = $(BUILDDIR)/lib
SRCDIR = fortran
MODDIR = $(SRCDIR)/modules

SRC = $(wildcard $(SRCDIR)/*.f90)
SRCFNAME = $(notdir $(SRC))
SRCBINS = $(SRCFNAME:%.f90=%)
SRCOBJECTS_O = $(addprefix $(INCLUDEDIR)/,$(addsuffix .obj,$(SRCBINS)))
SRCLIB = $(addsuffix .dll,$(SRCBINS))

MOD = $(wildcard $(MODDIR)/*.f90)
MODFNAME = $(notdir $(MOD))
MODBINS = $(MODFNAME:%.f90=%)
MODOBJECTS_O = $(addprefix $(INCLUDEDIR)/,$(addsuffix .obj,$(MODBINS)))
MODLIB = $(addprefix $(LIBDIR)/,$(addsuffix .dll,$(MODBINS)))

FFLAGS = -nologo -I$(INCLUDEDIR) -MD -Qparallel -Qmkl -Qopenmp -free -fpp -threads -heap-arrays- -integer-size:32 -dll -names:lowercase -nofixed -module:$(INCLUDEDIR)
# FFLAGS = -nologo -I$(INCLUDEDIR) -MD -Qparallel -fast -Qmkl -Qopenmp -free -fpp -threads -heap-arrays- -integer-size:32 -dll -names:lowercase -nofixed -module:$(INCLUDEDIR)

FLINKFLAGS = -link -nodefaultlib:libcmt

CFLAGS = -nologo -c -LD -MT -EHsc $(cdebug) -DMSDOS -D_CRT_SECURE_NO_WARNINGS

.ONESHELL:

.PHONY: clean sharedlib

sharedlib: $(SRCLIB)

$(INCLUDEDIR)/%.obj: $(MODDIR)/%.f90
	-mkdir -p $(INCLUDEDIR)
	ifort $(FFLAGS) -c $(MODDIR)/$(notdir $(@:.obj=.f90)) -o $@ $(FLINKFLAGS)

$(INCLUDEDIR)/cvputrec.obj: $(CODEVLIBDIR)/cvputrec.cpp
	-mkdir -p $(INCLUDEDIR)
	icl $(CFLAGS) $(CODEVLIBDIR)/cvputrec.cpp
	-mv -f cvputrec.obj $(INCLUDEDIR)

%.dll: $(SRCDIR)/%.f90 $(MODOBJECTS_O) $(INCLUDEDIR)/cvputrec.obj 
	-mkdir -p $(INCLUDEDIR)
	ifort $(FFLAGS) -o $@ $(SRCDIR)/$(notdir $(@:.dll=.f90)) $(MODOBJECTS_O) $(INCLUDEDIR)/cvputrec.obj $(FLINKFLAGS)
clean:
	rm -rf $(BUILDDIR)

