FEMDIR:=../..
include ../../Rules.make

EXEC = shyelab

INCLUDE:=-I.
INCLUDE+= -I../../fem3d/   # shypre requires param.h simul.h and netcdf.inc that are still in fem3d/
INCLUDE+= -I../../include/ # contains nosinf.h 

#--------------------------------------------------------------------------------------

ifneq ($(NETCDF),true)
  NETCDF_OBJ= \
	../../tools/elab/elab_nc_dummy.o
else
  LIB_NETCDF = -L$(NETCDFDIR)/lib -lnetcdf -L$(NETCDFFDIR)/lib -lnetcdff
  NETCDF_OBJ=\
	../../tools/elab/elab_nc.o \
	../../tools/2nc/netcdf.o \
	../../tools/2nc/netcdf_util.o 
endif


ifneq ($(PARALLEL_MPI),NONE)
  LIB_MPI = -L$(MPI_DIR)/lib -lmpi
  SHYMPI_OBJ=\
	../../src/mpi_node.o \
	../../fem3d/mod_shympi_node_internal.o \
	../../src/mpi_general.o 
else
  SHYMPI_OBJ=\
	../../src/mpi_dummy.o \
	../../src/mpi_general.o 
endif

#--------------------------------------------------------------------------------------

LFLAGS:= $(LIB_NETCDF) $(LIB_MPI)

#--------------------------------------------------------------------------------------

OBJ:=\
	../../src/basin.o \
	../../src/levels.o \
        $(SHYMPI_OBJ) \
	../../src/depth.o \
	../../src/depth_util.o \
	../../src/evgeom.o \
	../../tools/elab/elabtime.o \
	../../tools/elab/elab_off.o \
	../../src/geom_dynamic.o \
	../../src/clo.o \
	../../src/strings.o \
	../../tools/elab/elabutil.o \
	../../tools/elab/elab_gis.o \
	../../tools/elab/elab_expand.o \
	../../tools/elab/elab_regular.o \
	../../src/ioshy.o \
	../../src/mod_offline.o \
	../../tools/elab/shyelab_output.o \
	$(NETCDF_OBJ) \
	../../tools/elab/shyelab_util.o \
	../../src/round.o \
	../../src/regular.o \
	../../src/line_util.o \
	../../src/grd.o \
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
	../../src/zadaptutil.o \
	../../tools/elab/shyelab_extract.o \
	../../src/iso8601.o \
	../../src/dates.o \
	../../tools/elab/shyelab1.o \
	../../src/femfile.o \
	../../src/tsfile.o \
	../../src/ioflux.o \
	../../src/ioext.o \
	../../src/geometrical.o \
	../../src/fem_util.o \
	../../fem3d/subous.o \
	../../tools/2nc/output_util.o \
	../../tools/2nc/ousutil.o \
	../../src/mpi_buffer.o \
	../../src/sort_util.o \
	../../src/sigma.o \
	../../fem3d/subnos.o \
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
	../../src/to_linear.o \
	../../src/lagrange_read.o \
	../../tools/elab/lgrelab1.o \
	../../tools/elab/extelab1.o \
	../../tools/elab/flxelab1.o \
	../../src/util.o \
	../../src/heat_util.o \
	../../tools/elab/elab_resample.o \
	../../tools/elab/elab_check.o \
	../../tools/elab/tselab1.o \
	../../tools/elab/femelab1.o 

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
