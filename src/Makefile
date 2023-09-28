#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

FEMDIR = ..
include ../Rules.make

#--------------------	VARS	---------------------

EXEC = shyfem

#----------------------	INCLUDE	----------------------

INCLUDE = -I$(FEMDIR)/include -I$(PETSC_DIR)/include -I$(MPI_DIR)/include
INCLUDE = -I$(FEMDIR)/include 

#--------------------	GOTM	-----------------------
ifeq ($(GOTM),true)
  LIBF_GOTM      = $(LIBDIR)/libgotm.a
endif

#---------------------	PETSC	-----------------------
ifeq ($(SOLVER),PETSC)
  LIBG_PETSC = -L$(PETSC_DIR)/lib -lpetsc
endif

#---------------------	MPI	-----------------------
ifneq ($(PARALLEL_MPI),NONE)
  LIBG_MPI = -L$(MPI_DIR)/lib -lmpi
endif

#----------------------	LIBRARY	------------------------
LIBGS= $(LIBG_PETSC) $(LIBG_MPI)
LIBFS= $(LIBF_GOTM)

#----------------------------------------------------------------------------------------------------------------

OBJ = convert.o utility.o para.o iso8601.o time_util.o levels.o fil.o basin.o system_matrix.o \
      geom.o lnku.o tvd.o mpi_dummy.o mpi_ghost.o shyparts_dummy.o mpi_basin.o mpi_custom.o mpi_buffer.o mpi_time.o \
		mpi_general.o mpi_debug.o \
      links.o hydro_admin.o evgeom.o round.o mercury_dummy.o regular.o openmp_admin.o bnd_geom.o \
      meteo.o waves.o turbulence_util.o sinking.o nudging.o internal.o geom_dynamic.o depth.o \
      bnd_aux.o ts.o hydro_vel.o restart.o gotm_aux.o diff_aux.o bnd_dynamic.o area.o roughness.o diffusion.o hydro_baro.o \
      hydro_print.o dts.o femfile.o tsfile.o sigma.o intp.o defnames.o \
      version.o to_linear.o conz_util.o vert_intp.o nls.o bnd.o fem_util.o bnd_admin.o zadaptutil.o intp_fem_file.o \
      bnd_scalar.o projection.o coordinates.o subset.o tide_utils.o tidef.o layer_thickness.o bndo.o \
      conz_common.o bclfix_util.o nohyd.o outputd.o elems_dealing.o heat_util.o strings.o util.o meteo_forcing.o \
      heat_util2.o icemodel.o ice.o heat_admin.o topological.o heat_temp.o heat_mfsbulk.o heat_coare.o heat_gotm.o \
      heat_areg.o heat_gill.o heat_default.o chk_NaN.o meteo_admin.o heat_admin2.o sigma_admin.o \
      bit10.o ecological.o transforms.o tsfile_admin.o fluidmud.o mud_admin.o chezy.o grd.o line_admin.o \
      discharge_flux.o flux_util.o ioflux.o sedim_util.o waves_admin.o substress.o sedim_admin.o flxm.o flux.o geometrical.o volume.o ioext.o random.o depth_util.o \
      vec_util.o extra.o ets.o ioets.o ets_admin.o modls.o ioshy.o shy_util.o \
      ous_admin.o flux_box.o residual.o find.o tvd_admin.o check.o wetdry.o stability.o time_admin.o \
      model3d_util.o trace.o keps_util.o kepsilon.o eq_state.o mod_offline.o offline_data.o shy_turbulence.o rd_ts.o intp_tst.o \
      arrays.o grdutil.o closing.o system_spk.o nohydro.o diffusion_admin.o bcvel.o lagrange_data.o lagrange_sedim.o \
      lagrange_flux.o line_util.o lagrange_util.o random_gen.o int_to_cart.o lagrange_vertical.o lagrange_utils.o lagrange_connect.o \
      lagrange_read.o lagrange_inout.o lagrange_cont.o lagrange_larve.o lagrange_decay.o lagrange_init.o lagrange_util_tr.o \
      lagrange_track.o lagrange_dif.o ioeos.o lagrange_stk.o lagrange_main.o lagrange_back.o floodfill.o topological_admin.o \
      mpi_tripple.o explicit.o sort_util.o coo_matrix.o nudge.o spk_blas.o matrix_inv.o spk_matvec.o \
      spk_ilut.o spk_itaux.o spk_iters.o spk_coo.o sparskit_admin.o penta.o turbine.o hydrodynamic.o submir.o initialize.o \
      bnd_routines.o concentration_omp.o histogram.o concentration.o baroclinic.o def_para.o \
      iostr.o constants.o f_vol.o mkdistance.o conz_admin.o dates.o water_ren_time.o toxical.o custom_admin.o \
      befor_after.o bfm_dummy.o bfm.o test_zeta.o clo.o lnkmpi.o poisson.o ww3_dummy.o simple_sedi.o 

      #mod_petsc_global.o mod_petsc_system.o simsys_petsc.o
#----------------------------------------------------------------------------------------------------------------

%.o: %.f 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

%.o: %.f90 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

$(EXEC):$(MODULES) $(OBJ) $(LIBFS) shyfem.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) -o $@ $^ $(LIBGS) 
	cp $(EXEC) $(BINDIR)
	ar r $(LIBDIR)/libshy.a $(OBJ)
clean:
	-rm -f $(MODDIR)/*.mod *.o
	-rm -f $(EXEC)
	-rm -f $(BINDIR)/$(EXEC)
	-rm -f $(LIBDIR)/libshy.a 
cleanall: clean

#----------------------------------------------------------------------------------------------------------------
