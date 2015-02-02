FC = gfortran
FFLAGS=-cpp -ffixed-line-length-none -fno-align-commons -g
NCINC=-I/usr/include
NCLIBS=-lnetcdff -lnetcdf
PROGS=gen-coef convert-sst-from-nc write-swr convert-albedo-ascii-to-nc \
      convert-forest-fraction-ascii-to-nc convert-orog-ascii-to-dat

all: $(PROGS)

gen-coef: gen-coef.f90
	$(FC) $(FFLAGS) -I../src -o $@ $^

convert-sst-from-nc: convert-sst-from-nc.f90
	$(FC) $(FFLAGS) $(NCINC) -I../src -o $@ $^ $(NCLIBS)

write-swr: write-swr.f
	$(FC) $(FFLAGS) -o $@ $^

convert-albedo-ascii-to-nc: convert-albedo-ascii-to-nc.f90
	$(FC) $(FFLAGS) $(NCINC) -I../src -o $@ $^ $(NCLIBS)

convert-forest-fraction-ascii-to-nc: convert-forest-fraction-ascii-to-nc.f90
	$(FC) $(FFLAGS) $(NCINC) -I../src -o $@ $^ $(NCLIBS)

convert-orog-ascii-to-dat: convert-orog-ascii-to-dat.f90
	$(FC) $(FFLAGS) $(NCINC) -I../src -o $@ $^ $(NCLIBS)

clean:
	rm -f $(PROGS)
