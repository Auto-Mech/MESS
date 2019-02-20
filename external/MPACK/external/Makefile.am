SUBDIRS = dummy

if !SYSTEM_GMP
SUBDIRS += gmp
endif

if !SYSTEM_MPFR
SUBDIRS += mpfr
endif

if !SYSTEM_MPC
SUBDIRS += mpc
endif

if !SYSTEM_QD
SUBDIRS += qd
endif

all:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i || exit 1; \
	done

clean:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i clean || exit 1; \
	done
	rm -rf i
