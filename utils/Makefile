FC = gfortran
FFLAGS=-cpp -ffixed-line-length-none -fno-align-commons -g
NCINC=-I/usr/include
NCLIBS=-lnetcdff -lnetcdf
PROGS=tonc convert-sst-ascii-to-nc convert-sst-nc-to-ascii tograds

all: $(PROGS)

tonc: tonc.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

convert-sst-ascii-to-nc: convert-sst-ascii-to-nc.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

convert-sst-nc-to-ascii: convert-sst-nc-to-ascii.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

tograds: tograds.f
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f $(PROGS)
