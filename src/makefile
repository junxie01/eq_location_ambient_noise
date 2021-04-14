FC=gfortran
FFLAG=-ffixed-line-length-none
subjects1=locat.o sort.o cal_dist.o
subjects2=locat_onev.o sort.o cal_dist.o
subjects3=locat_vofcrust1.0.o sort.o cal_dist.o
subjects4=locat_chi.o sort.o cal_dist.o
subjects5=bootstrapping_locate.o sort.o cal_dist.o init_random_seed.o
subjects6=locate_with_all4station.o sort.o cal_dist.o 
all:../bin/locat ../bin/locat_onev ../bin/locat_vofcrust1.0 ../bin/locat_chi ../bin/locat_bootstrap ../bin/locate_with_all4station
%.o:%.f
	$(FC) -c $(FFLAG) $^
%.o:%.f90
	$(FC) -c $(FFLAG) $^
../bin/locat:$(subjects1)
	$(FC) $(subjects1) -o $@ $(FFLAG)
../bin/locat_onev:$(subjects2)
	$(FC) $(subjects2) -o $@ $(FFLAG)
../bin/locat_vofcrust1.0:$(subjects3)
	$(FC) $(subjects3) -o $@ $(FFLAG)
../bin/locat_chi:$(subjects4)
	$(FC) $(subjects4) -o $@ $(FFLAG)
../bin/locat_bootstrap:$(subjects5)
	$(FC) $^ -o $@ $(FFLAG)
../bin/locate_with_all4station:$(subjects6)
	$(FC) $^ -o $@ $(FFLAG)
clean:
	rm $(subjects1) $(subjects2) $(subjects3) $(subjects4)
