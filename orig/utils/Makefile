FC = gfortran
FFLAGS=-cpp -ffixed-line-length-none -fno-align-commons -g
NCINC=-I/usr/include
NCLIBS=-lnetcdff -lnetcdf
PROGS=inputfiles tonc convert-sst-ascii-to-nc convert-sst-ascii-to-dat \
      convert-sst-nc-to-ascii convert-orog-ascii-to-nc tograds

all: $(PROGS)

inputfiles: inputfiles.f
	$(FC) $(FFLAGS) -I../src -o $@ $^

tonc: tonc.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

convert-sst-ascii-to-nc: convert-sst-ascii-to-nc.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

convert-sst-ascii-to-dat: convert-sst-ascii-to-dat.f
	$(FC) $(FFLAGS) -I../src -o $@ $^

convert-sst-nc-to-ascii: convert-sst-nc-to-ascii.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

convert-orog-ascii-to-nc: convert-orog-ascii-to-nc.f90
	$(FC) $(FFLAGS) $(NCINC) -o $@ $^ $(NCLIBS)

tograds: tograds.f
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f $(PROGS)
