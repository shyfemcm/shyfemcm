SHELL:= /bin/bash
FC:=mpif90
INCLUDE:=-I.
INCLUDE+= -I../../fem3d/   # shypre requires param.h simul.h and netcdf.inc that are still in fem3d/
INCLUDE+= -I../../include/ # contains nosinf.h 

FLAGS:= $(INCLUDE) -cpp 

LFLAGS:= -L$(NETCDF_C_HOME) -lnetcdf -L$(NETCDF_FORTRAN_HOME) -lnetcdff

OBJ:=\
	../../src/basin.o \
	../../src/levels.o \
	../../fem3d/mod_shympi_node.o \
	../../fem3d/mod_shympi_node_internal.o \
	../../fem3d/shympi_general.o \
	../../src/depth.o \
	../../src/depth_util.o \
	../../tools/2nc/netcdf.o \
	../../src/geom_dynamic.o \
	../../src/evgeom.o \
	../../tools/pre/subclo.o \
	../../src/strings.o \
	../../tools/elab/elabutil.o \
	../../tools/elab/shyelab_util.o \
	../../src/ioshy.o \
	../../fem3d/mod_offline.o \
	../../tools/elab/elabtime.o \
	../../tools/elab/elab_off.o \
	../../tools/elab/shyelab_output.o \
	../../tools/2nc/netcdf_util.o \
	../../tools/elab/elab_nc.o \
	../../tools/elab/elab_gis.o \
	../../tools/elab/elab_expand.o \
	../../tools/elab/elab_regular.o \
	../../fem3d/subround.o \
	../../src/regular.o \
	../../fem3d/line_util.o \
	../../tools/pre/grd.o \
	../../tools/elab/shyelab_area.o \
	../../tools/elab/shyelab_average.o \
	../../tools/elab/shyelab_nodes.o \
	../../src/coordinates.o \
	../../src/projection.o \
	../../tools/elab/shyelab_proj.o \
	../../src/chk_NaN.o \
	../../src/vec_util.o \
	../../src/shy_util.o \
	../../fem3d/shyutil.o \
	../../fem3d/subshyutil2.o \
	../../fem3d/zadaptutil.o \
	../../tools/elab/shyelab_extract.o \
	../../src/iso8601.o \
	../../fem3d/subdates.o \
	../../tools/elab/shyelab1.o \
	../../src/femfile.o \
	../../src/tsfile.o \
	../../src/ioflux.o \
	../../src/ioext.o \
	../../src/geometrical.o \
	../../src/ioous.o \
	../../src/fem_util.o \
	../../tools/2nc/output_util.o \
	../../tools/2nc/ousutil.o \
	../../fem3d/submpi_buffer.o \
	../../src/sort_util.o \
	../../src/sigma.o \
	../../src/ionos.o \
	../../tools/2nc/nosutil.o \
	../../src/dts.o \
	../../fem3d/unitutils.o \
	../../fem3d/subdum.o \
	../../src/utility.o \
	../../src/para.o \
	../../src/nls.o \
	../../src/util/subnsa.o \
	../../src/def_para.o \
	../../src/defnames.o \
	../../src/fil.o \
	../../src/util/apnfile.o \
	../../src/convert.o \
	../../src/version.o \
	../../fem3d/sub2linear.o \
	../../fem3d/lagrange_read.o \
	../../tools/elab/lgrelab1.o \
	../../tools/elab/extelab1.o \
	../../tools/elab/flxelab1.o \
	../../fem3d/subutil.o \
	../../src/heat_util.o \
	../../tools/elab/elab_resample.o \
	../../tools/elab/elab_check.o \
	../../tools/elab/tselab1.o \
	../../tools/elab/femelab1.o \
	../../tools/elab/shyelab.o 

%.o:%.f90
	${FC} $(FLAGS) -c $< -o $@ 

shyelab: $(OBJ)
	$(FC) $(LFLAGS) -o shyelab $(OBJ)

clean:
	rm ../../*/*/*.o
	rm ../../*/*.o
	rm ../../*/*/*.mod
	rm ../../*/*.mod
	rm shyelab
