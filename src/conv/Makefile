BIN=${LIPSIA_BIN}/.

#SUBDIRS = vnifti vdicom vprtheader vold2new vcopygeo vaddgeo vscalevals
SUBDIRS = $(shell ls -d */)

all::
	for i in $(SUBDIRS:/=) ;\
	do \
		(cd $$i ; make && cp $$i ${BIN}); \
	done


clean::
	for i in $(SUBDIRS/=) ;\
	do \
		(cd $$i ; make clean); \
	done
