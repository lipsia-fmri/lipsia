BIN=${LIPSIA_BIN}/.

#SUBDIRS = vcompcor vdenoise vdespike vpreprocess vslicetime vsrad vrepairstart vdetrend
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
