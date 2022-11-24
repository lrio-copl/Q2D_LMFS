CODEVLIBDIR = C:/CODEV114/umr/lib

BUILDDIR = fortran/build

BINDIR = $(BUILDDIR)/bin
INCLUDEDIR = $(BUILDDIR)/lib
LIBDIR = $(BUILDDIR)/lib
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

FFLAGS = -nologo -I$(INCLUDEDIR) -MD -Qparallel -Qmkl -Qopenmp -free -fpp -threads -heap-arrays- -Qipo- -integer-size:32 -dll -names:lowercase -nofixed -module:$(INCLUDEDIR)
# FFLAGS = -nologo -I$(INCLUDEDIR) -MD -Qparallel -fast -Qmkl -Qopenmp -free -fpp -threads -heap-arrays- -integer-size:32 -dll -names:lowercase -nofixed -module:$(INCLUDEDIR)

FLINKFLAGS = -link -nodefaultlib:libcmt

# FFLAGS2 = -nologo -I$(INCLUDEDIR) -MT -Qparallel -Qmkl -Qopenmp -free -fpp -threads -heap-arrays- -Qipo- -integer-size:32 -dll -names:lowercase -nofixed -module:$(INCLUDEDIR) 

MKL_LINKS = mkl_intel_c.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib 
# -nodefaultlib:libcmt

CFLAGS = -nologo -c -LD -MT -EHsc $(cdebug) -DMSDOS -D_CRT_SECURE_NO_WARNINGS

.ONESHELL:

.PHONY: clean sharedlib

sharedlib: $(SRCLIB) mkl.dll

mkl.dll:
	@echo Build MKL custom dll
	mkdir -p $(LIBDIR)
	MKL_LIB="mkl_intel_c_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib" 
	cd "$(MKLROOT)/tools/builder";nmake ia32 export=$(subst \c\,C:\,$(subst /,\\,$(shell pwd)/$(SRCDIR)/user.def)) name=mkl;cp mkl.dll $(subst \c\,c:\,$(subst /,\\,$(shell pwd)));cp mkl.lib $(subst \c\,c:\,$(subst /,\\,$(shell pwd)));cd -
	cp "$(MKLROOT)/../../compiler/latest/windows/compiler/lib/ia32/libiomp5md.lib" $(LIBDIR)/libiomp5md.lib
	cp "$(MKLROOT)/../../compiler/latest/windows/compiler/lib/ia32/libiomp5md.lib" .
	cp "$(MKLROOT)/../../compiler/latest/windows/compiler/lib/ia32/libiomp5md_db.dll" libiomp5md.dll

$(INCLUDEDIR)/%.obj: $(MODDIR)/%.f90
	-mkdir -p $(INCLUDEDIR)
	ifort $(FFLAGS) -c $(MODDIR)/$(notdir $(@:.obj=.f90)) -o $@ $(FLINKFLAGS)

$(INCLUDEDIR)/cvputrec.obj: $(CODEVLIBDIR)/cvputrec.cpp
	-mkdir -p $(INCLUDEDIR)
	icl $(CFLAGS) $(CODEVLIBDIR)/cvputrec.cpp
	-mv -f cvputrec.obj $(INCLUDEDIR)

%.dll: $(SRCDIR)/%.f90 $(MODOBJECTS_O) $(INCLUDEDIR)/cvputrec.obj mkl.dll
	-mkdir -p $(INCLUDEDIR)
	ifort $(FFLAGS) -o $@ $(SRCDIR)/$(notdir $(@:.dll=.f90)) $(MODOBJECTS_O) $(INCLUDEDIR)/cvputrec.obj $(FLINKFLAGS)
clean:
	rm -rf $(BUILDDIR)

