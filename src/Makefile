# $LIPSIA_DEV/src/Makefile
# /bk2015
#
# CFLAGS, LDFLAGS, CPPFLAGS are set using "$LIPSIA_DEV/lipsia-setup.sh"
# so sourcing this script is mandatory
#

SUBDIRS  = $(shell ls -d */)
VIA_LIBS = lib_viaio lib_via lib_niftiio

# VIA_LIBS are a global prequisite, so we make the dependency explicit
# here and below, to enable parallel builds like "make -j 4"
.PHONY: subdirs $(VIA_LIBS) $(SUBDIRS) install clean binclean distclean

ifeq ($(and $(strip $(LIPSIA_LIB)),$(strip $(LIPSIA_INCLUDE)),$(strip $(LIPSIA_INST)),$(strip $(LIPSIA_BIN))),)
  $(error "Please at least run "source ../lipsia-setup.sh" before !" )
endif

subdirs: $(SUBDIRS)

$(VIA_LIBS):
	$(MAKE) -C $@

$(SUBDIRS): $(VIA_LIBS)
	$(MAKE) -C $@

clean:
	$(foreach D, $(SUBDIRS), $(MAKE) -C $D clean;)
	find . -name \*.o -delete

binclean:
	rm -rf ../lib/* ../bin/*

distclean: binclean clean

install: subdirs
	echo installing ...
	cp ../bin/* $(LIPSIA_INST)/bin
	cp ../lib/* $(LIPSIA_INST)/lib

###
