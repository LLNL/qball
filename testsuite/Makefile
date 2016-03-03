SUBDIRS = 01_hydrogen 02_nitrogen

.PHONY: check all $(SUBDIRS)

all:
check: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ check

