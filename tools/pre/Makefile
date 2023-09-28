FEMDIR:=../..
include ../../Rules.make

EXEC = shypre

INCLUDE:=-I.
INCLUDE+= -I../../fem3d/   # shypre requires param.h and simul.h that are still in fem3d/

#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------

LFLAGS:= 

#--------------------------------------------------------------------------------------

OBJ:=\
	../../tools/pre/rosen.o \
	../../src/basin.o \
	../../src/fem_util.o \
	../../tools/pre/cuthill_mckee.o \
	../../src/geom_dynamic.o \
	../../tools/pre/grd.o \
	../../tools/pre/grdutil.o \
	../../src/defnames.o \
	../../fem3d/subdum.o \
	../../src/utility.o \
	../../src/para.o \
	../../src/nls.o \
	../../src/strings.o \
	../../src/sort_util.o \
	../../src/util/subnsa.o \
	../../src/def_para.o \
	../../src/convert.o \
	../../src/fil.o \
	../../src/util/apnfile.o \
	../../src/version.o \
	../../tools/pre/subclo.o \
	../../fem3d/subsort.o 
%.o:%.f
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.f90
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.c
	${CC} $(CFLAGS) $(INCLUDE) -c $< -o $@

$(EXEC):$(MODULES) $(OBJ) $(LIBFS) $(EXEC).o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBGS) 
	cp $(EXEC) ../../bin/	

clean:
	-rm $(FEMDIR)/*/*.o 
	-rm $(FEMDIR)/*/*/*.o
	-rm $(FEMDIR)/*/*/*/*.o
	-rm $(FEMDIR)/lib/mod/*
	-rm $(EXEC)

