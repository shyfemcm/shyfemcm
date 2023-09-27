SHELL:= /bin/bash
FC:=ifort
INCLUDE:=../../fem3d/  # shypre requires param.h and simul.h that are still in fem3d/
FLAGS:= -I. -I$(INCLUDE) -cpp

OBJ:= \
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
	../../fem3d/subsort.o \
	../../tools/pre/shypre.o 

%.o:%.f90
	${FC} $(FLAGS) -c $< -o $@ 

shypre: $(OBJ)
	$(FC) -o shypre $(OBJ)

clean:
	rm *.mod
	rm ../../tools/pre/rosen.o 
	rm ../../src/basin.o 
	rm ../../src/fem_util.o 
	rm ../../tools/pre/cuthill_mckee.o 
	rm ../../src/geom_dynamic.o 
	rm ../../tools/pre/grd.o 
	rm ../../tools/pre/grdutil.o 
	rm ../../src/defnames.o 
	rm ../../fem3d/subdum.o 
	rm ../../src/utility.o 
	rm ../../src/para.o 
	rm ../../src/nls.o 
	rm ../../src/strings.o 
	rm ../../src/sort_util.o 
	rm ../../src/util/subnsa.o 
	rm ../../src/def_para.o 
	rm ../../src/convert.o 
	rm ../../src/fil.o 
	rm ../../src/util/apnfile.o 
	rm ../../src/version.o 
	rm ../../tools/pre/subclo.o 
	rm ../../fem3d/subsort.o 
	rm ../../tools/pre/shypre.o 
	rm shypre 
