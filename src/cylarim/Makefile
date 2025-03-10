BIN=~/bin

SUBDIRS = vmetric vrim vcylarim vcylarim_stats vcylarim_getmask vcylarim_plot vcylarim_seed vcylarim_cover
#SUBDIRS = $(shell ls -d */)

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
