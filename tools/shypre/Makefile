
FEMDIR:=../..
include $(FEMDIR)/Rules.make

EXEC = shypre

INCLUDE = -I$(FEMDIR)/include/

LIBLUT  = -L$(DIRLIB) -l$(LIBUTIL)
LIBFUT  = $(DIRLIB)/lib$(LIBUTIL).a

LIBLSU  = -L$(DIRLIB) -l$(LIBSHYUTIL)
LIBFSU  = $(DIRLIB)/lib$(LIBSHYUTIL).a

LIBL = $(LIBLPS) $(LIBLUT) $(LIBLSU)
LIBF = $(LIBFPS) $(LIBFUT) $(LIBFSU)

#----------------------------------------------------------------------------

LFLAGS:= 

#----------------------------------------------------------------------------

OBJ:=\
	../../tools/pre/rosen.o \
	../../tools/pre/cuthill_mckee.o \
	../../src/basin.o \
	../../src/fem_util.o \
	../../src/geom_dynamic.o \
	../../src/grd.o \
	../../src/grdutil.o \
	../../src/defnames.o \
	../../src/utility.o \
	../../src/para.o \
	../../src/nls.o \
	../../src/strings.o \
	../../src/sort_util.o \
	../../src/def_para.o \
	../../src/convert.o \
	../../src/fil.o \
	../../src/version.o \
	../../src/clo.o \
	../../fem3d/subsort.o 

#	../../src/util/apnfile.o \
#	../../src/util/subnsa.o \

#----------------------------------------------------------------------------

$(EXEC): $(OBJ) $(LIBF) $(EXEC).o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBL) 
	cp $(EXEC) ../../bin/	

clean:
	-rm $(FEMDIR)/*/*.o 
	-rm $(FEMDIR)/*/*/*.o
	-rm $(FEMDIR)/*/*/*/*.o
	-rm $(FEMDIR)/lib/mod/*
	-rm $(EXEC)

cleanall: clean

#----------------------------------------------------------------------------

%.o:%.f
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.f90
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.c
	${CC} $(CFLAGS) $(INCLUDE) -c $< -o $@

#----------------------------------------------------------------------------

