builddir = build
outputdir = output
liseqsrcdir = src/liseq

forsrcdir = src/fortran
formoddir = $(forsrcdir)/modules
forsrcfiles = $(wildcard $(forsrcdir)/*.f90) $(wildcard $(formoddir)/*.f90)

codevfigeps = $(wildcard $(builddir)/fig/*_plt.eps)
codevfigepsname = $(notdir $(codevfigeps))
codevfigpdf = $(codevfigepsname:%_plt.eps=%.pdf)
completecodevpdf = $(addprefix $(outputdir)/fig/, $(codevfigpdf))

codevliseq = $(wildcard $(liseqsrcdir)/*.liseq) $(wildcard $(liseqsrcdir)/*/*.liseq)
codevseq = $(codevliseq:$(liseqsrcdir)/%.liseq=$(builddir)/%.seq)


V = 1
ACTUAL_CC = codev-remote 
CC_0 = @ $(ACTUAL_CC)
CC_1 = $(ACTUAL_CC) -v
CC = $(CC_$(V))

.PHONY: clean all sync

.ONESHELL:

all : $(builddir)/fortran_stamp $(builddir)/codev_stamp $(completecodevpdf)

clean: 
	-@ rm -r $(outputdir) $(builddir)
	$(CC) -c

$(codevseq): $(builddir)/%.seq: $(liseqsrcdir)/%.liseq 
	-@ mkdir -p $(builddir)/$(subst $(liseqsrcdir)/,,$(dir $<))
	-@ liseq $< -o $@

$(builddir)/fortran_stamp : $(forsrcfiles)
	-@ mkdir -p $(builddir)
	-@ mkdir -p $(outputdir)
	$(CC) -k --config_path codev_remote_build.yaml 
	$(CC) -d --config_path codev_remote_build.yaml 
	$(CC) -c --config_path codev_remote_build.yaml 
	-@ rsync -aP $(outputdir)/*.dll $(builddir)/ -q
	-@ touch $(builddir)/fortran_stamp

$(builddir)/codev_stamp : $(codevliseq)  $(codevseq) $(forsrcfiles)
	-@ mkdir -p $(builddir)
	-@ mkdir -p $(outputdir)
	-@ cp media/*.seq $(builddir)
	$(CC)
	# $(CC) -c --config_path codev_remote_build.yaml 
	-@ touch $(builddir)/codev_stamp

$(completecodevpdf) : %.pdf : %_plt.eps 
	-@ codev-style-vie $< > $(outputdir)/fig/temp.eps
	-@ pstopdf $(outputdir)/fig/temp.eps -o $@ -p
	-@ mogrify -transparent-color white -density 600 -format png $@
	-@ rm $(outputdir)/fig/temp.eps

sync : $(wildcard */*.seq) $(wildcard *.seq) $(forsrcfiles) $(codevseq)
	-@ mkdir -p $(builddir)
	-@ mkdir -p $(outputdir)
	-@ rsync -aP $(builddir)/bin/. ./ -q
	$(CC) -p 
	$(CC) -d 
	-@ rsync -aP output/ $(builddir) -q
	touch $(builddir)/codev_stamp
	-@ rm -r $(codevliseq:%.liseq=%.seq)
	# $(MAKE)
