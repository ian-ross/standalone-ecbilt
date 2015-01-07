FC = gfortran
FFLAGS=-cpp -ffixed-line-length-none -fno-align-commons -g

SRCS=ecbilt.f global.f initial.f atmmodel.f atmphys.f atmdyn.f atmdiag.f \
     landmodel.f fluxmodel.f oceanfixed.f error.f nag.f root.f

INCS=comatfor.h comatm.h comcoup.h comcouphelp.h comdiag.h comdyn.h \
     comglobal.h comice.h comlake.h comcouplake.h comland.h comocean.h \
     comocfix.h comoutlocal.h comphys.h openatinpfiles.h openatoutfiles.h \
     openocfixfiles.h openatstartfiles.h

OBJS=$(addprefix src/, $(SRCS:.f=.o))
INCALL=$(addprefix src/, $(INCS))

ecbilt :$(OBJS)
	$(FC) $(FFLAGS) -o ecbilt $(OBJS)

src/%.o: src/%.f $(INCALL)
	$(FC) $< $(FFLAGS) -c -o $@

clean:
	rm -f $(OBJS)
	rm -f ecbilt
