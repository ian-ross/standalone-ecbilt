NETCDFINC=-I/usr/local/netcdf-cxx-4.1.3/include
NETCDFLIB=-L/usr/local/netcdf-cxx-4.1.3/lib -lnetcdf -lnetcdff
FC=gfortran
FINCS=-I./include $(NETCDFINC)
FFLAGS=-Wall -cpp -ffixed-line-length-none -fno-align-commons $(FINCS)
LD=gfortran
LDFLAGS=$(NETCDFLIB)

SRCDIR=./src
OBJDIR=./obj

BASESRCS=emic.f coupling0.f ecbilt0.f initial0.f atmphys0.f atmdyn0.f \
         atmoutp0.f atmdiag0.f nag.f sst_pert.f error0.f
SRCS=$(addprefix $(SRCDIR)/, $(BASESRCS))
OBJS=$(addprefix $(OBJDIR)/, $(BASESRCS:.f=.o))

standalone-ecbilt: $(OBJDIR) $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) ${LIBS}

obj:
	if [ ! -d obj ]; then mkdir obj; fi

clean:
	rm -f $(OBJDIR)/*.o
	rm -f standalone-ecbilt

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	$(FC) $< $(FFLAGS) -c -o $@
