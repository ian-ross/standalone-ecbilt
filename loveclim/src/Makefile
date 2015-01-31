NETCDFINCLUDE=-I/usr/local/netcdf-cxx-4.1.3/include
NETCDFLIB=-L/usr/local/netcdf-cxx-4.1.3/lib -lnetcdf -lnetcdff
UDUNITSINCLUDE=-I/usr/include
UDUNITSLIB=-L/usr/lib -ludunits

FC=gfortran
FFLAGS=-ffixed-line-length-none -fno-align-commons -Jobj ${NETCDFINCLUDE}
LD=gfortran
LIBS=${NETCDFLIB} ${UDUNITSLIB}

SRCmain=main.f90
SRCATMOS=ecbilt.f90 initial.f90 atmphys.f90 atmdyn.f90 atmoutp.f90 \
         atmdiag.f90 nag.f coupling.f90 error.f90
SRCOCEAN=oceanfixed.f90 root.f
SRCLAND=landmodel.f90

INCALL=comatfor.h comatm.h comcoup.h comcouphelp.h comdiag.h comdyn.h \
       comemic.h comland.h comoutlocal.h comphys.h openatinpfiles.h \
       openatoutfiles.h comocfix.h

SRCALL=$(SRCmain) $(SRCATMOS) $(SRCLAND) $(SRCOCEAN)
OBJALL1=$(SRCALL:.f90=.o)
OBJALL=$(addprefix obj/, $(OBJALL1:.f=.o))

ecbilt: obj $(OBJALL)
	$(LD) $(LDFLAGS) -o  $@ $(OBJALL) ${LIBS}
	chmod g+rwx $@

obj:
	if [ ! -d obj ]; then mkdir obj; fi

$(OBJmain): $(INCALL)

$(OBJATMOS): $(INCALL)

$(OBJLAND): $(INCALL)

obj/%.o: %.f90
	$(FC) $< $(FFLAGS) -c -o $@

obj/%.o: %.f
	$(FC) $< $(FFLAGS) -c -o $@

clean:
	rm -fr ecbilt obj
