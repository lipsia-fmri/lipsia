
BIN=${LIPSIA_BIN}/.

#SUBDIRS = vconjunction vonesample_ttest vpaired_ttest vpaired_ranktest vtwosample_ttest vtwosample_ranktest vtwosample_ttest_tfce vonesample_ttest_tfce vtwosample_wilcox
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
