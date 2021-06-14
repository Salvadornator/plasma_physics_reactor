# to create libinterposg95.a (does not need link to lapack routines, includes dpgbtrf_s.f90)
PROG          = interposg95
# to create libinterposifort.a
PROG_ifort     = interposifort
# to create libinterpospgi.a
PROG_pgi     = interpospgi
# to create libinterpospgi.a
PROG_pgf95     = interpospgf95
# to create libinterpospgi.a
PROG_gfortran  = interposgfortran
# to create libinterposg95_wolapack.a (needs link to lapack routines)
PROG_wolapack = interposg95_wolapack
# to create libinterposifort_wolapack.a (needs link to lapack routines)
PROG_ifortwolapack = interposifort_wolapack
# to create libinterpospgi_wolapack.a (needs link to lapack routines)
PROG_pgiwolapack = interpospgi_wolapack
PROG_pgf95wolapack = interpospgf95_wolapack
PROG_gfortranwolapack = interposgfortran_wolapack
# to create interpostest_f90 program from interpostest.f90
PROG_interpostest_f90 = interpostest_f90
PROG_interpostest_f90_pgi = interpostest_f90_pgi
PROG_interpostest_eqdsk = interpostest_f90_eqdsk
# to compile matlab interpos routine: either from within matlab with mex interpos.f90 after creating interpos.f90
# create interpos.f90 for matlab including already lapack routines
PROG_interpos_for_matlab = interpos_for_matlab
PROG_interpos_for_matlab32bits = interpos_for_matlab32bits
# create interpos.f90 for matlab and needs link to wolapack
PROG_interpos_for_matlab_wolapack = interpos_for_matlab_wolapack
# mex directly matlab routine with own lapack routine, no link needed
PROG_interpos = interpos_mex
# mex directly matlab routine with own lapack routine, no link needed
PROG_interpos_wolapack = interpos_mex_wolapack
# interpos_source.f90 for CHEASE or other standalone versions (avoid having a library linked)
PROG_interpos_source = interpos_source

SRCS_LIB =	cbsplgen.f90 cbsplgenp.f90 cbfitbnd.f90 cbfitper.f90 \
		splipera.f90 splibnda.f90 intlinear.f90 intquadratic.f90 
MODS =	prec_rkind.f90 interpos_module.f90
SRCS_topmatlab = interpos_matlab_top.f90
SRCS_topmatlab32bits = interpos_matlab_top32bits.f90
SRCS_testf90 = interpostest_profile.f90
SRCS_testf90 = interpostest3.f90
SRCS_testf90_eqdsk = interpostest_2Deqdsk.f90 indexx.f90

SRCS_lap = dpgbtrf_s.f90

OBJS_LIB = $(SRCS_LIB:.f90=.o)

LIBS_nolap =	
LIBS_lap =	-L${MKL}   -lmkl_lapack -lmkl -lguide -lpthread
F90_95 =	g95
F90_ifort =	ifort
F90_pgi =	pgf90
F90_pgf95 =	pgf95
F90_gfortran =	gfortran
FC_95 =	-fPIC -fno-omit-frame-pointer
FC_ifort =	-O2 -r8 -fPIC
FC_pgi	=	 -r8 -fPIC -Mnosecond_underscore
FC_pgf95	=	 -r8 -fPIC -Mnosecond_underscore
FC_gfortran = -g -O0 -Wall -ffree-line-length-none

all1:	$(PROG) $(PROG_ifort) $(PROG_pgi) $(PROG_pgf95) $(PROG_wolapack) $(PROG_ifortwolapack) $(PROG_pgiwolapack) $(PROG_pgf95wolapack) $(PROG_interpos_for_matlab) $(PROG_interpos_source) $(PROG_gfortran) $(PROG_gfortranwolapack) 
$(PROG):
	$(F90_95) -c $(FC_95) $(MODS)
	$(F90_95) -c $(FC_95) $(SRCS_LIB) $(SRCS_lap)
	ar -rv lib$@.a  $(MODS:.f90=.o) $(OBJS_LIB) $(SRCS_lap:.f90=.o)
	mv $(MODS:.f90=.mod) lib$@.a g95_libs
	echo "now can compile with: $(F90_95) -o file_a_out $(FC_95) file.f90 -L$(PWD) -l$@"

$(PROG_ifort):
	$(F90_ifort) -c $(FC_ifort) $(MODS)
	$(F90_ifort) -c $(FC_ifort) $(SRCS_LIB) $(SRCS_lap)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB) $(SRCS_lap:.f90=.o)
	mv $(MODS:.f90=.mod) lib$@.a ifort_libs
	echo "now can compile with: $(F90_ifort) -o file_a_out $(FC_ifort) file.f90 -L$(PWD) -l$@"

$(PROG_pgi):
	$(F90_pgi) -c $(FC_pgi) $(MODS)
	$(F90_pgi) -c $(FC_pgi) $(SRCS_LIB) $(SRCS_lap)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB) $(SRCS_lap:.f90=.o)
	mv $(MODS:.f90=.mod) lib$@.a pgi_libs
	echo "now can compile with: $(F90_pgi) -o file_a_out $(FC_pgi) file.f90 -L$(PWD) -l$@"

$(PROG_pgf95):
	$(F90_pgf95) -c $(FC_pgf95) $(MODS)
	$(F90_pgf95) -c $(FC_pgf95) $(SRCS_LIB) $(SRCS_lap)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB) $(SRCS_lap:.f90=.o)
	mv $(MODS:.f90=.mod) lib$@.a pgf95_libs
	echo "now can compile with: $(F90_pgf95) -o file_a_out $(FC_pgf95) file.f90 -L$(PWD)-l$@"

$(PROG_gfortran):
	$(F90_gfortran) -c $(FC_gfortran) $(MODS)
	$(F90_gfortran) -c $(FC_gfortran) $(SRCS_LIB) $(SRCS_lap)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB) $(SRCS_lap:.f90=.o)
	mv $(MODS:.f90=.mod) lib$@.a gfortran_libs
	echo "now can compile with: $(F90_gfortran) -o file_a_out $(FC_gfortran) file.f90 -L$(PWD)-l$@"

$(PROG_wolapack):
	$(F90_95) -c $(FC_95) $(MODS)
	$(F90_95) -c $(FC_95) $(SRCS_LIB)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB)
	mv lib$@.a g95_libs
	echo "now can compile with: $(F90_95) -o file_a_out $(FC_95) file.f90 -L$(PWD) $(LIBS_lap) -l$@"

$(PROG_ifortwolapack):
	$(F90_ifort) -c $(FC_ifort) $(MODS)
	$(F90_ifort) -c $(FC_ifort) $(SRCS_LIB)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB)
	mv lib$@.a ifort_libs
	echo "now can compile with: $(F90_ifort) -o file_a_out $(FC_ifort) file.f90 -L$(PWD) $(LIBS_lap) -l$@"

$(PROG_pgiwolapack):
	$(F90_pgi) -c $(FC_pgi) $(MODS)
	$(F90_pgi) -c $(FC_pgi) $(SRCS_LIB)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB)
	mv lib$@.a pgi_libs
	echo "now can compile with: $(F90_pgi) -o file_a_out $(FC_pgi) file.f90 -L$(PWD) $(LIBS_lap) -l$@"

$(PROG_pgf95wolapack):
	$(F90_pgf95) -c $(FC_pgf95) $(MODS)
	$(F90_pgf95) -c $(FC_pgf95) $(SRCS_LIB)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB)
	mv $(MODS:.f90=.mod) lib$@.a pgf95_libs
	echo "now can compile with: $(F90_pgf95) -o file_a_out $(FC_pgf95) file.f90 -L$(PWD) $(LIBS_lap) -l$@"

$(PROG_gfortranwolapack):
	$(F90_gfortran) -c $(FC_gfortran) $(MODS)
	$(F90_gfortran) -c $(FC_gfortran) $(SRCS_LIB)
	ar -rv lib$@.a $(MODS:.f90=.o) $(OBJS_LIB)
	mv $(MODS:.f90=.mod) lib$@.a gfortran_libs
	echo "now can compile with: $(F90_gfortran) -o file_a_out $(FC_gfortran) file.f90 -L$(PWD) $(LIBS_lap) -l$@"

$(PROG_interpostest_f90):
	$(F90_95) -o $@ $(FC_95) $(SRCS_testf90) -I$(PWD)/g95_libs -L$(PWD)/g95_libs -l$(PROG)

$(PROG_interpostest_f90_pgi): $(PROG_pgi)
	$(F90_pgi) -o $@ $(FC_pgi) $(SRCS_testf90) -L$(PWD) -l$(PROG_pgi)

$(PROG_interpostest_eqdsk):
	$(F90_gfortran) -o $@ $(FC_gfortran) $(SRCS_testf90_eqdsk)  -I$(PWD)/gfortran_libs -L$(PWD)/gfortran_libs -linterposgfortran
#	$(F90_95) -o $@ $(FC_95) $(SRCS_testf90_eqdsk)  -I$(PWD)/g95_libs -L$(PWD)/g95_libs -l$(PROG)
#	$(F90_pgi) -o $@ $(FC_pgi) $(SRCS_testf90_eqdsk) -I$(PWD)/pgi_libs -L$(PWD)/pgi_libs -l$(PROG_pgi)

$(PROG_interpos_for_matlab): 
	cat $(MODS) $(SRCS_topmatlab) $(SRCS_LIB) $(SRCS_lap) > ../matlab/interpos.f90
	echo "within matlab >> mex interpos.f90"

$(PROG_interpos_for_matlab32bits): 
	cat $(MODS) $(SRCS_topmatlab32bits) $(SRCS_LIB) $(SRCS_lap) > ../matlab/interpos32bits.f90
	echo "after cp interpos32bits.f90 interpos.f90, within matlab >> mex interpos.f90"

$(PROG_interpos_for_matlab_wolapack): 
	cat $(MODS) $(SRCS_topmatlab) $(SRCS_LIB) > ../matlab/interpos.f90
	echo "within matlab >> mex FLIBS='-Wl,-rpath-link,/usr/local/matlab-7.7.0/bin/glnxa64 -L/usr/local/matlab-7.7.0/bin/glnxa64 -lmx -lmex -lmat -lm -lmwlapack' interpos.f90"

$(PROG_interpos):	$(PROG_interpos_for_matlab)
	g95 -c  -I/usr/local/matlab-7.7.0/extern/include -I/usr/local/matlab-7.7.0/simulink/include -fexceptions -fPIC -fno-omit-frame-pointer  -DMX_COMPAT_32 -O3  ../matlab/interpos.f90
	g95 -O3 -pthread -shared -Wl,--version-script,/usr/local/matlab-7.7.0/extern/lib/glnxa64/fexport.map -Wl,--no-undefined -o  "../matlab/interpos.mexa64"  interpos.o  "/usr/local/matlab-7.7.0/extern/lib/glnxa64/version4.o"  -Wl,-rpath-link,/usr/local/matlab-7.7.0/bin/glnxa64 -L/usr/local/matlab-7.7.0/bin/glnxa64 -lmx -lmex -lmat -lm
	rm interpos.o

$(PROG_interpos_wolapack):	$(PROG_interpos_for_matlab_wolapack)
	g95 -c  -I/usr/local/matlab-7.7.0/extern/include -I/usr/local/matlab-7.7.0/simulink/include -fexceptions -fPIC -fno-omit-frame-pointer  -DMX_COMPAT_32 -O3  ../matlab/interpos.f90
	g95 -O3 -pthread -shared -Wl,--version-script,/usr/local/matlab-7.7.0/extern/lib/glnxa64/fexport.map -Wl,--no-undefined -o  "../matlab/interpos.mexa64"  interpos.o  "/usr/local/matlab-7.7.0/extern/lib/glnxa64/version4.o"  -Wl,-rpath-link,/usr/local/matlab-7.7.0/bin/glnxa64 -L/usr/local/matlab-7.7.0/bin/glnxa64 -lmx -lmex -lmat -lm -lmwlapack
	rm interpos.o

$(PROG_interpos_source): 
	cat $(MODS) $(SRCS_LIB) > interpos_source.f90
	@echo "you can now copy interpos_source.f90 to your directory"

clean:	
	rm *.o *.a *.mod

source:	
	cat $(MODS) $(SRCS_LIB) > ../interpos_source.f90

