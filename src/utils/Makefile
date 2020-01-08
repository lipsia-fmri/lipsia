BIN=${LIPSIA_BIN}/.

#SUBDIRS = vattrcopy vattredit vave vcheckdata vsubsequence vstats vsymmetric vtimestep vshowmosaic vflip3d vmerge
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
