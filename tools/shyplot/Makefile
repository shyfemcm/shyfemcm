FEMDIR:=../..
include ../../Rules.make

EXEC = shyplot

#--------------------------------------------------------------------------------------

INCLUDE:=-I.
INCLUDE+= -I../../fem3d/   # shypre requires param.h simul.h and netcdf.inc that are still in fem3d/
INCLUDE+= -I../../include/ # contains nosinf.h 
INCLUDE+= -I../../tools/postscript/ # contains psgraph.h

#--------------------------------------------------------------------------------------


ifneq ($(NETCDF),true)
  NETCDF_OBJ:= \
	../../tools/elab/elab_nc_dummy.o 
else
  LIB_NETCDF = -L$(NETCDFDIR)/lib -lnetcdf -L$(NETCDFFDIR)/lib -lnetcdff
  NETCDF_OBJ=\
	../../tools/2nc/netcdf.o \
	../../tools/2nc/netcdf_util.o \
	../../tools/elab/elab_nc.o 
endif

ifneq ($(PARALLEL_MPI),NONE)
  LIB_MPI = -L$(MPI_DIR)/lib -lmpi
  SHYMPI_OBJ=\
	../../src/mpi_node.o \
	../../src/mpi_node_internal.o \
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
	../../src/geom_dynamic.o \
	../../src/evgeom.o \
	../../src/clo.o \
	../../src/strings.o \
	../../tools/elab/elabutil.o \
	../../tools/elab/shyelab_util.o \
	../../src/ioshy.o \
	../../src/mod_offline.o \
	../../tools/elab/elabtime.o \
	../../tools/elab/elab_off.o \
	../../tools/elab/shyelab_output.o \
	$(NETCDF_OBJ) \
	../../tools/elab/elab_gis.o \
	../../tools/elab/elab_expand.o \
	../../tools/elab/elab_regular.o \
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
	../../fem3d/subous.o \
	../../src/fem_util.o \
	../../tools/2nc/output_util.o \
	../../tools/2nc/ousutil.o \
	../../src/mpi_buffer.o \
	../../src/sort_util.o \
	../../src/sigma.o \
	../../fem3d/subnos.o \
	../../tools/2nc/nosutil.o \
	../../src/dts.o \
	../../fem3d/unitutils.o \
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
	../../tools/elab/femelab1.o \
	../../tools/shyplot/mod_hydro_plot.o \
	../../src/geom.o \
	../../src/hydro_admin.o \
	../../tools/postscript/general.o \
	../../tools/postscript/gustd.o \
	../../src/hydro_vel.o \
	../../tools/shyplot/colortb.o \
	../../tools/shyplot/supcust.o \
	../../tools/postscript/colorutil.o \
	../../tools/postscript/psgraph.o \
	../../tools/shyplot/util/col.o \
	../../tools/shyplot/util/supcol.o \
	../../tools/shyplot/suptim.o \
	../../src/intp.o \
	../../src/time_util.o \
	../../src/vert_intp.o \
	../../tools/shyplot/util/suputi.o \
	../../src/intp_fem_file.o \
	../../tools/shyplot/supint.o \
	../../tools/shyplot/supbas.o \
	../../tools/postscript/psgraphf.o \
	../../tools/shyplot/supano.o \
	../../src/lnku.o \
	../../src/bnd.o \
	../../src/bnd_geom.o \
	../../src/links.o \
	../../src/bnd_admin.o \
	../../src/hydro_print.o \
	../../src/topological.o \
	../../fem3d/color.o \
	../../tools/shyplot/plotutil.o \
	../../src/grdutil.o \
	../../src/ioeos.o \
	../../fem3d/submask.o \
	../../tools/shyplot/supdep.o \
	../../fem3d/mod_connect.o \
	../../tools/shyplot/supver.o \
	../../tools/shyplot/suplin1.o \
	../../tools/shyplot/supout.o \
	../../tools/shyplot/suplin.o \
	../../tools/shyplot/supiso.o \
	../../tools/shyplot/supsim.o \
	../../tools/shyplot/suplag.o

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
