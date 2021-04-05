#Remover Tool
RM = /bin/rm -rf

#Copy Tool
CP = /bin/cp -a

#Make Directory tool
MKDIR = mkdir

#Source Directory
SRCDIR = src

#Library Name
LFILE = liblrt.a

#Top-level rule.
all: lrtlib

lrtlib:
	$(RM) $(LFILE)
	@ for f in $(SRCDIR); do (cd $$f && echo In $$f && make all); done
	ln -s $(SRCDIR)/$(LFILE) $(LFILE)

install: lrtlib
ifndef LRTPATH
	@echo
	@echo -----------------------------------------------------------
	@echo You need to define LRTPATH in order to install the program.
	@echo -----------------------------------------------------------
	@echo
else 
ifeq ($(PWD),$(LRTPATH)) 
	@echo
	@echo -----------------------------------------------------------
	@echo LRTPATH points to the current folder. 
	@echo No need to proceed with the installation.
	@echo -----------------------------------------------------------
	@echo
else
	@ [ -d $(LRTPATH) ] || $(MKDIR) $(LRTPATH)
	$(CP) Filters specs $(LRTPATH)
	$(CP) $(SRCDIR)/$(LFILE) $(LRTPATH) 
	@ [ -d $(LRTPATH)/src ] || $(MKDIR) $(LRTPATH)/src
endif	
endif


uninstall:
ifndef LRTPATH
	@echo
	@echo -----------------------------------------------------------
	@echo You need to define LRTPATH in order to uninstall the program.
	@echo -----------------------------------------------------------
	@echo
else 
ifeq ($(PWD),$(LRTPATH)) 
	@echo
	@echo -----------------------------------------------------------
	@echo LRTPATH points to the current folder. 
	@echo Cannot uninstall.
	@echo -----------------------------------------------------------
	@echo
else
	$(RM) -r $(LRTPATH)
endif
endif


distclean: clean
	$(RM) $(LFILE) $(SRCDIR)/$(LFILE)
clean: 
	@ for s in $(SRCDIR); do (cd $$s && echo In $$s && make clean ); done

