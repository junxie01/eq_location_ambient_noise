FC=gfortran
FFLAG=-ffixed-line-length-none
subjects=bootstrapping_locate.o sort.o cal_dist.o init_random_seed.o
all:bootstrapping_locate
%.o:%.f
	$(FC) -c $(FFLAG) $^
%.o:%.f90
	$(FC) -c $(FFLAG) $^
bootstrapping_locate:$(subjects)
	$(FC) $(subjects) -o $@ $(FFLAG)
clean:
	rm $(subjects)
