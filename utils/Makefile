FC = gfortran
FFLAGS=-cpp -ffixed-line-length-none -fno-align-commons -g
NCINC=-I/usr/include
NCLIBS=-lnetcdff -lnetcdf
PROGS=tonc convert-sst tograds

all: $(PROGS)

tonc: tonc.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

convert-sst: convert-sst.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

tograds: tograds.f
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f $(PROGS)
