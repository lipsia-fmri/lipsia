BIN=${LIPSIA_BIN}/.

MFLAGS =

#SUBDIRS = vbinarize vcontrast vdelsmall vgauss3d vhemi vlabel3d vquickmorph3d vselbig vsmooth3d vdist3d
SUBDIRS = $(shell ls -d */)

all::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS:/=) ;\
	do \
		(cd $$i ; make $(MFLAGS) ; cp $$i ${BIN}); \
	done


clean::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS:/=) ;\
	do \
		(cd $$i ; make clean); \
	done
