FC=gfortran
FFLAG=-ffixed-line-length-none
subjects=locat.o sort.o cal_dist.o
all:locat
install:
	cp locat ../bin
%.o:%.f
	$(FC) $^ -c $(FFLAG)
%.o:%.f90
	$(FC) $^ -c $(FFLAG)
locat:$(subjects)
	$(FC) $(subjects) -o $@ $(FFLAG)
clean:
	rm $(subjects)
