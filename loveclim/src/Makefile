NETCDFINCLUDE=-I/usr/local/netcdf-cxx-4.1.3/include
NETCDFLIB=-L/usr/local/netcdf-cxx-4.1.3/lib -lnetcdf -lnetcdff
UDUNITSINCLUDE=-I/usr/include
UDUNITSLIB=-L/usr/lib -ludunits

FC=gfortran
FFLAGS=-ffixed-line-length-none -fno-align-commons ${NETCDFINCLUDE}
LD=gfortran
LIBS=${NETCDFLIB} ${UDUNITSLIB}

SRCmain=emic.f
SRCATMOS=ecbilt0.f initial0.f atmphys0.f atmdyn0.f atmoutp0.f atmdiag0.f nag.f \
         coupling0.f error0.f
SRCVEGET=veget.f veget_outp.f veget_sub.f veget_wr.f
SRCLAND=landmodel0.f
SRCSIO=clio.f start.f flucor.f scale.f slopes.f scadew.f scadns.f scali.f \
       uve.f uvi.f barot.f uvb0et.f uvbfet.f uvm.f informe.f defcst.f \
       defgrid.f geogra.f redforc.f correct.f cortab.f conti3d.f etat.f \
       vdiffu.f slopez.f alph2dc.f alphdkc.f alphdec.f raccord.f inforun.f\
       savrunb.f redrunb.f savtab.f redtab.f foroutp.f rahmflx.f ec_oc2co.f \
       ec_co2oc.f sloptab.f
SRCISO=isoslope.f isoadvec.f isodiffu.f
SRCTKE=diftur.f engtur.f
SRCICE=staocb.f outave.f ocean_outp.f outnet.f streamfunc.f streaminit.f \
       forcat.f shine.f scater.f gather.f thersf.f acrlbq.f fontbc.f \
       ocesla.f forcng.f pdecli.f icdyna.f icdadv.f dynami.f zmaxt.f advx.f \
       advy.f diffus.f adv.f
SRCOCEAN=$(SRCSIO) $(SRCICE) $(SRCTKE) $(SRCISO)

INCALL=comatfor.h comatm.h comcoup.h comcouphelp.h comdiag.h comdyn.h \
       comemic.h comland.h comoutlocal.h comphys.h openatinpfiles.h \
       openatoutfiles.h veget.h

INCSIO=bloc.com comclio.h datadc.com ice.com para.com thermo.com vareq.com \
       bloc0.com densit.com isoslope.com para0.com trace.com varno.com \
       const.com dynami.com moment.com reper.com type.com

SRCEMIC=$(SRCmain) $(SRCATMOS) $(SRCLAND) $(SRCOCEAN) $(SRCVEGET)
OBJEMIC=$(addprefix obj/, $(SRCEMIC:.f=.o))

emic.x: obj $(OBJEMIC)
	$(LD) $(LDFLAGS) -o  $@ $(OBJEMIC) ${LIBS}
	chmod g+rwx $@

obj:
	if [ ! -d obj ]; then mkdir obj; fi

$(OBJmain): $(INCALL)

$(OBJATMOS): $(INCALL)

$(OBJLAND): $(INCALL)

$(OBJOCEAN): $(INCALL) $(INCSIO)

obj/%.o: %.f
	$(FC) $< $(FFLAGS) -c -o $@

clean:
	rm -f emic.x $(OBJEMIC)
