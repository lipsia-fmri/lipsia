BIN=${LIPSIA_BIN}

#SUBDIRS = vncm vccm vpaired_ccm vcorr vecm vspectralecm veta2 vpowerspectrum vReHo vspectralcluster vtc
SUBDIRS = $(shell ls -d */)

all::
	for i in $(SUBDIRS:/=) ;\
	do \
		(cd $$i ; make && cp $$i ${BIN}); \
	done


clean::
	for i in $(SUBDIRS:/=) ;\
	do \
		(cd $$i ; make clean); \
	done
