#----------------------------------------------
# makefile for Chek_MPI
# Giorgio Amati, CASPUR 2000-2012
# Giorgio Amati, CINECA 2013-20??
#----------------------------------------------
# notes:
#
# 	1) .....
#
#----------------------------------------------
# targets:
#
# 
# flags
#
#----------------------------------------------
#
DATE = $(shell date +%d%m%y)
#
# default
#----------------------------------------------
FINE = "********** Fine Compilazione *********"
#----------------------------------------------
#
#----------------------------------------------
# compiler section
#----------------------------------------------
# default (NVIDIAcompiler)
 CC = mpicc
 FC = mpifort
 FIX = -DPGI
 FOPT = -fast -Mcontiguous -Mnodepchk
 OPT = -fast -Mcontiguous -Mnodepchk
 mpiP = -L /leonardo/home/userinternal/gamati01/LOCAL_SW/mpiP/lib -lmpiP -lm -lbfd
#
# laptop
ifdef GNU
 CC = mpicc
 FC = mpif90
 FOPT = -ftree-vectorize -Ofast
 OPT  = -ftree-vectorize -Ofast
 mpiP = -L ~/LOCAL_SW/bin/mpiP/lib/ -lmpiP -lm -lbfd
endif
#
# Galileo
ifdef INTEL
 CC = mpiicc
# FC = mpiifort
 FC = mpiifx
# FOPT = -O3 -xCORE-AVX512 -mtune=skylake-avx512  -assume contiguous_pointer
# OPT  = -O3 -xCORE-AVX512 -mtune=skylake-avx512  -assume contiguous_pointer
 FOPT = -O3 -assume contiguous_pointer
 OPT  = -O3 -assume contiguous_pointer
 mpiP =
endif
#
#----------------------------------------------
# version section
#----------------------------------------------
#
# default
VER = step3
#
FIX=-DSTEP3
GPU= -acc -Minfo=acc -DOPENACC -DFAST
VER=step3
#
ifdef STEP7
 FIX= -DSTEP7
 GPU= -acc -Minfo=acc -DOPENACC -DFAST
 VER=step7
endif
#
ifdef STEP9
 FIX= -DSTEP9
 GPU= -acc -Minfo=acc -DOPENACC -DFAST
 OBJ1= bcond.comm.step9.o \
       col_MC_masked.o
 VER=step9
endif
#
ifdef STEP10
 FIX= -DSTEP10
 GPU= -acc -Minfo=acc -DOPENACC -DFAST
 OBJ1= bcond.comm.step10.o \
       do_somethingGPU_masked.o
 VER=step10
endif
#
#----------------------------------------------
# Other stuff
#----------------------------------------------
#
ifdef GPUENABLE
 FIX+= -DGPUENABLE -DOPENACC
 FOPT+=$(GPU) 
endif
#
ifdef MPIP
 FIX+= -DMPIP
 LIB =$(mpiP) 
endif
#----------------------------------------------
# Makefile section
#----------------------------------------------
MOD  = get_mem_mod.o \
       precision.o \

OBJ  = storage.o \
       mod_tools.o \
       input.o \
       alloca.o \
       bgk3d_HIT.o \
       col_MC.o \
       finalize.o \
       diagno.o \
       diagnostic.o \
       dissipation.o \
       do_somethingGPU.o \
       hencol.o \
       init.o \
       init_old.o \
       boundaries.o \
       bcond.comm.step0.o \
       bcond.comm.step1.o \
       bcond.comm.step2.o \
       bcond.comm.step3.o \
       bcond.comm.step4.o \
       bcond.comm.step5.o \
       bcond.comm.step6.o \
       bcond.comm.step7.o \
       setup_MPI.o \
       probe_global.o \
       prof_i.o \
       prof_j.o \
       prof_k.o \
       profile.o \
       git_info.o \
       make_info.o \
       outdat.o \
       vtk_xy.binary.o \
       vtk_xz.binary.o \
       vtk_yz.binary.o \
       vtk_3d.binary.o \
#
#----------------------------------------------
#
all: $(MOD) $(OBJ) $(OBJ1)
	$(FC) $(FOPT) $(OBJ) $(OBJ1) $(FIX) $(TRY) $(LIB) -o ../RUN/bgk3d_HIT.$(VER).x 
	rm make_info.F90
	rm git_info.F90
	@echo $(FINE)


#----------------------------------------------

%.o %.mod: $(INC) %.f90
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(TRY) -c $<

%.o: $(INC) %.F90
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(TRY) -c $<

%.o: $(INC) %.f
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(TRY) c $<

%.o: $(INC) %.c
	$(CC) $(COPT)                          -c $<

%.p: %.f90
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(TRY) $<

#----------------------------------------------
make_info.F90:
	echo "    subroutine make_info"                      >  make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FIX)    ` \"  >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FC)     ` \"  >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FOPT)   ` \"  >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(PREPROC)` \"  >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(GPU)    ` \"  >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(TRY)    ` \"  >> make_info.F90
	echo "    end subroutine make_info"                  >> make_info.F90

git_info.F90:
	echo "    subroutine git_info"                              >> git_info.F90
	echo "    write( 6,*) \""INFO: `git log  | head  -n 1` \"   >> git_info.F90
	echo "    write(16,*) \""INFO: `git log  | head  -n 1` \"   >> git_info.F90
	echo "    end subroutine git_info"                          >> git_info.F90
#
#----------------------------------------------
clean:
	rm -f *.o;
	rm -f *.mod;
	rm -f core;
#----------------------------------------------
info:  $(INC)
	@echo "Modules          =  "$(MOD);
	@echo "Objects          =  "$(OBJ);
	@echo "Compiler         =  "$(FC);
	@echo "Compiler flags   =  "$(FOPT);
	@echo "GPU      flags   =  "$(GPU);
	@echo "Preproc  flags   =  "$(PREPROC);
	@echo "Fix      flags   =  "$(FIX);
	@echo "Try      flags   =  "$(TRY);
	@echo "Double   flags   =  "$(DBL);
	@echo "Architecture     =  "$(HOSTTYPE);
	@echo "Version          =  "$(VER);
	@echo "Libraries        =  "$(LIB);

#----------------------------------------------
tar:
	make clean
	rm -rf ../RUN/*.vtk
	rm -rf ../RUN/*.bin
	rm -rf ../RUN/*.h5
	tar -cvf ../../archive.bgk_HIT.$(DATE).tar ../../MPI_GPUAWARE/*
#----------------------------------------------
help:
	@echo "mpi (default)"
#----------------------------------------------

