
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

FEMDIR = ../..
include $(FEMDIR)/Rules.make
FEMBIN = $(FEMDIR)/bin

#--------------------	EXES	---------------------

EXES = shyfem

#----------------------	INCLUDE	----------------------

INCLUDE = -I$(FEMDIR)/include
ifeq ($(SOLVER),PETSC)
  INCLUDE += -I$(PETSC_DIR)/include
endif

#--------------------	LIBS	-----------------------

LIBUTIL         = util
LIBSHYUTIL      = shyutil
LIBSHYMPI       = shympi
LIBSHYUTILMPI   = shyutilmpi
LIBSPARSKIT     = sparskit

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBLMP  = -L$(LIBDIR) -l$(LIBSHYMPI)
LIBFMP  = $(LIBDIR)/lib$(LIBSHYMPI).a

LIBLUM  = -L$(LIBDIR) -l$(LIBSHYUTILMPI)
LIBFUM  = $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBLSP  = -L$(LIBDIR) -l$(LIBSPARSKIT)
LIBFSP  = $(LIBDIR)/lib$(LIBSPARSKIT).a

LIBF = $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT) $(LIBFSP)
LIBL = $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIBLGOTM) $(LIBG_METIS)
LIBL += $(LIBLGOTM) $(LIBG_METIS) $(LIBG_PETSC) $(LIBLSP)

#--------------------	GOTM	-----------------------

ifeq ($(GOTM),true)
  LIBLGOTM      = -L $(LIBDIR) -lgotm
  LIBFGOTM      = $(LIBDIR)/libgotm.a
endif

#---------------------	SOLVER	-----------------------

SOLVER_OBJ = simsys_lp.o 

ifeq ($(SOLVER),SPARSKIT)
   SOLVER_OBJ = \
	simsys_spkaux.o sparskit_admin.o coo_util.o \
	spk_ilut.o spk_iters.o spk_itaux.o spk_coo.o \
	spk_itpack.o spk_matvec.o spk_blas.o matrix_inv.o
   SOLVER_OBJ = \
	simsys_spkaux.o sparskit_admin.o coo_util.o \
	matrix_inv.o
endif

ifeq ($(SOLVER),PETSC)
   LIBG_PETSC = -L$(PETSC_DIR)/lib -lpetsc
   SOLVER_OBJ		 = spk_blas.o matrix_inv.o \
	amgx_c_wrapper_dummy.o mod_petsc_global.o \
	mod_petsc_system.o simsys_petsc.o
endif

ifeq ($(PARTS),METIS)
  LIBG_METIS = -L$(METISDIR)/lib -lmetis
endif

#------------------------------------------------------------------------

LAGRANGE_OBJ = \
      lagrange_data.o lagrange_sedim.o \
      lagrange_flux.o lagrange_util.o \
      lagrange_vertical.o lagrange_utils.o lagrange_connect.o \
      lagrange_inout.o lagrange_cont.o lagrange_larve.o \
      lagrange_decay.o lagrange_init.o lagrange_util_tr.o \
      lagrange_track.o lagrange_dif.o lagrange_stk.o \
      lagrange_main.o lagrange_back.o \

HEAT_OBJ = \
      meteo.o meteo_forcing.o meteo_admin.o \
      icemodel.o ice.o \
      waves.o waves_admin.o ww3_dummy.o \
      heat_util2.o heat_admin.o heat_admin2.o \
      heat_temp.o heat_mfsbulk.o heat_coare.o heat_gotm.o \
      heat_areg.o heat_gill.o heat_default.o \

ECO_OBJ = \
      ecological.o bfm_dummy.o bfm.o mercury_dummy.o \

SEDIMENT_OBJ = \
      sedim_util.o simple_sedi.o sedim_admin.o \

MPI_EXTRA_OBJ = \
      mpi_tripple.o mpi_assert.o lnkmpi.o

OBJ = system_matrix.o tvd.o \
      hydro_admin.o openmp_admin.o \
      turbulence_util.o sinking.o nudging.o internal.o \
      ts.o hydro_vel.o restart.o gotm_aux.o diff_aux.o \
      area.o roughness.o diffusion.o hydro_baro.o \
      hydro_print.o conz_util.o \
      bnd_scalar.o subset.o tide_utils.o tidef.o \
      layer_thickness.o bndo.o conz_common.o \
      bclfix_util.o nohyd.o outputd.o elems_dealing.o sigma_admin.o \
      transforms.o tsfile_admin.o fluidmud.o mud_admin.o chezy.o \
      line_admin.o discharge_flux.o flux_util.o \
      flxm.o modls.o \
      ous_admin.o residual.o tvd_admin.o \
      wetdry.o stability.o time_admin.o \
      model3d_util.o trace.o keps_util.o kepsilon.o eq_state.o \
      offline_data.o rd_ts.o intp_tst.o \
      arrays.o closing.o diffusion_admin.o bcvel.o \
      floodfill.o topological_admin.o coo_matrix.o nudge.o \
      turbine.o initialize.o \
      bnd_routines.o concentration_omp.o histogram.o \
      concentration.o baroclinic.o def_para.o \
      f_vol.o mkdistance.o conz_admin.o \
      water_ren_time.o atoxi.o toxical.o \
      befor_after.o iostr.o constants.o \
      test_zeta.o \
	$(SOLVER_OBJ) \
	$(HEAT_OBJ) \
	substress.o \
	$(SEDIMENT_OBJ) \
	$(LAGRANGE_OBJ) \
	$(ECO_OBJ) \
	$(MPI_EXTRA_OBJ) \
	turbulence.o nohydro.o explicit.o hydrodynamic.o \
	custom_admin.o poisson.o \
      flux.o extra.o boxes.o check.o \

OBSOLETE = \
	bnd_aux.o bnd_dynamic.o tsfile.o bnd_admin.o \
	mod_offline.o \

#------------------------------------------------------------------------

default: fem

fem: $(EXES) links nothing

$(EXES): $(MODULES) $(OBJ) $(LIBF) shyfem.o
	 $(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

depend:
	@$(FEMBIN)/mkdp.pl -moddir=$(MODDIR) *.f90

nothing:
	@true

#------------------------------------------------------------------------

%.o: %.f 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

%.o: %.f90 
	$(F77) $(FFLAGS) $(INCLUDE) -c $< 

clean:
	-rm -f *.o
	-rm -f $(EXES)
	-rm -f $(BINDIR)/$(EXES)

cleanall: clean
	-rm -f $(MODDIR)/*.mod
	-rm -f $(LIBDIR)/libshy.a 

#update:
#	cp $(EXES) /home/georg/work/shyfem_repo/shyfem-mpi/bin
#	cp $(EXES) /home/georg/work/shyfem_repo/shyfem-mpi/fem3d
#	cp $(EXES) /home/georg/work/shyfem_repo/shyfem/bin
#	cp $(EXES) /home/georg/work/shyfem_repo/shyfem/fem3d

#------------------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../lib/mod/arrays.mod: arrays.o
../../lib/mod/befor_after.mod: befor_after.o
../../lib/mod/box.mod: boxes.o
../../lib/mod/box_arrays.mod: boxes.o
../../lib/mod/check_unit.mod: check.o
../../lib/mod/chezy.mod: chezy.o
../../lib/mod/close.mod: closing.o
../../lib/mod/connectivity.mod: lagrange_connect.o
../../lib/mod/extra.mod: extra.o
../../lib/mod/flux.mod: flux.o
../../lib/mod/histo.mod: histogram.o
../../lib/mod/ice_global.mod: icemodel.o
../../lib/mod/ice_params.mod: icemodel.o
../../lib/mod/ice_variables.mod: icemodel.o
../../lib/mod/lgr_sedim_module.mod: lagrange_sedim.o
../../lib/mod/lines.mod: line_admin.o
../../lib/mod/meteo_forcing_module.mod: meteo_forcing.o
../../lib/mod/mod_area.mod: area.o
../../lib/mod/mod_bclfix.mod: bclfix_util.o
../../lib/mod/mod_bfm.mod: bfm.o
../../lib/mod/mod_bndo.mod: bndo.o
../../lib/mod/mod_bstress.mod: substress.o
../../lib/mod/mod_conz.mod: conz_common.o
../../lib/mod/mod_diff_aux.mod: diff_aux.o
../../lib/mod/mod_diff_visc_fric.mod: diffusion.o
../../lib/mod/mod_fluidmud.mod: fluidmud.o
../../lib/mod/mod_gotm_aux.mod: gotm_aux.o
../../lib/mod/mod_hydro.mod: hydro_admin.o
../../lib/mod/mod_hydro_baro.mod: hydro_baro.o
../../lib/mod/mod_hydro_print.mod: hydro_print.o
../../lib/mod/mod_hydro_vel.mod: hydro_vel.o
../../lib/mod/mod_internal.mod: internal.o
../../lib/mod/mod_keps.mod: turbulence_util.o
../../lib/mod/mod_lagrange.mod: lagrange_data.o
../../lib/mod/mod_layer_thickness.mod: layer_thickness.o
../../lib/mod/mod_meteo.mod: meteo.o
../../lib/mod/mod_nohyd.mod: nohyd.o
../../lib/mod/mod_nudge.mod: nudge.o
../../lib/mod/mod_nudging.mod: nudging.o
../../lib/mod/mod_parwaves.mod: waves_admin.o
../../lib/mod/mod_petsc_global.mod: mod_petsc_global.o
../../lib/mod/mod_petsc_system.mod: mod_petsc_system.o
../../lib/mod/mod_renewal_time.mod: water_ren_time.o
../../lib/mod/mod_restart.mod: restart.o
../../lib/mod/mod_roughness.mod: roughness.o
../../lib/mod/mod_sediment.mod: sedim_admin.o
../../lib/mod/mod_sediment_para.mod: sedim_admin.o
../../lib/mod/mod_sedtrans05.mod: sedim_util.o
../../lib/mod/mod_shyfem.mod: shyfem.o
../../lib/mod/mod_sinking.mod: sinking.o
../../lib/mod/mod_subset.mod: subset.o
../../lib/mod/mod_system.mod: system_matrix.o
../../lib/mod/mod_system_global.mod: simsys_petsc.o
../../lib/mod/mod_system_global.mod: simsys_spkaux.o
../../lib/mod/mod_system_global.mod: system_spk.o
../../lib/mod/mod_system_interface.mod: system_matrix.o
../../lib/mod/mod_test_zeta.mod: test_zeta.o
../../lib/mod/mod_ts.mod: ts.o
../../lib/mod/mod_tvd.mod: tvd.o
../../lib/mod/mod_waves.mod: waves.o
../../lib/mod/poisson.mod: poisson.o
../../lib/mod/semi_implicit.mod: iostr.o
../../lib/mod/shyice_model.mod: ice.o
../../lib/mod/shympi_tripple.mod: mpi_tripple.o
../../lib/mod/simple_sediments.mod: simple_sedi.o
../../lib/mod/stab.mod: stability.o
../../lib/mod/subqfx.mod: heat_admin.o
../../lib/mod/subww3.mod: ww3_dummy.o
../../lib/mod/tide.mod: tide_utils.o
amgx_c_wrapper_dummy.o: ../../lib/mod/petscdm.mod \
		../../lib/mod/petscdmlabel.mod \
		../../lib/mod/petscmat.mod \
		../../lib/mod/petscvec.mod 
baroclinic.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/pkonst.mod \
		../../lib/mod/shympi.mod 
bcvel.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_bclfix.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/shympi.mod 
bfm.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_sinking.mod \
		../../lib/mod/mod_ts.mod 
bnd_routines.o: ../../lib/mod/basin.mod \
		../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod ../../lib/mod/mkonst.mod \
		../../lib/mod/mod_bnd_aux.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/pkonst.mod ../../lib/mod/shympi.mod 
bnd_scalar.o: ../../lib/mod/intp_fem_file.mod \
		../../lib/mod/iso8601.mod ../../lib/mod/shympi.mod 
bndo.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/shympi.mod 
bndo_admin.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_bndo.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/shympi.mod 
boxes.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_debug.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod \
		../../lib/mod/simul.mod 
check.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/femtime.mod ../../lib/mod/levels.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_nohyd.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod 
chezy.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_fluidmud.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/nls.mod ../../lib/mod/pkonst.mod \
		../../lib/mod/shympi.mod 
closing.o: ../../lib/mod/arrays.mod ../../lib/mod/basin.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/shympi.mod 
concentration.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod ../../lib/mod/mkonst.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_aux.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod \
		../../lib/mod/shympi_debug.mod 
concentration_omp.o: ../../lib/mod/basin.mod \
		../../lib/mod/evgeom.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_aux.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_subset.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod 
constants.o: ../../lib/mod/basin.mod ../../lib/mod/befor_after.mod \
		../../lib/mod/mkonst.mod ../../lib/mod/mod_bnd.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/pkonst.mod 
conz_admin.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/para.mod \
		../../lib/mod/shympi.mod 
conz_common.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_restart.mod \
		../../lib/mod/shympi.mod 
conz_util.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_ts.mod 
coo_matrix.o: ../../lib/mod/mod_system.mod
custom_admin.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/femtime.mod \
		../../lib/mod/iso8601.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_fluidmud.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_ts.mod \
		../../lib/mod/mod_waves.mod \
		../../lib/mod/shympi.mod 
def_para.o: ../../lib/mod/para.mod
diffusion_admin.o: ../../lib/mod/basin.mod \
		../../lib/mod/evgeom.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_aux.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/shympi.mod 
discharge_flux.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_debug.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_vel.mod 
elems_dealing.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/shympi.mod 
explicit.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/pkonst.mod \
		../../lib/mod/shympi.mod \
		../../lib/mod/shympi_debug.mod \
		../../lib/mod/shympi_tripple.mod 
extra.o: ../../lib/mod/basin.mod ../../lib/mod/befor_after.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_sediment.mod \
		../../lib/mod/mod_ts.mod \
		../../lib/mod/mod_waves.mod ../../lib/mod/nls.mod \
		../../lib/mod/shympi.mod ../../lib/mod/simul.mod 
f_vol.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_layer_thickness.mod 
floodfill.o: ../../lib/mod/basin.mod ../../lib/mod/mod_geom.mod
flux.o: ../../lib/mod/basin.mod ../../lib/mod/befor_after.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_sediment.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/nls.mod \
		../../lib/mod/shympi.mod ../../lib/mod/simul.mod 
flux_util.o: ../../lib/mod/levels.mod ../../lib/mod/mod_geom.mod
flxm.o: ../../lib/mod/shympi.mod ../../lib/mod/simul.mod
gotm_aux.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_restart.mod 
heat_admin2.o: ../../lib/mod/basin.mod \
		../../lib/mod/heat_const.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_debug.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_ts.mod 
heat_areg.o: ../../lib/mod/heat_const.mod
heat_coare.o: ../../lib/mod/heat_const.mod
heat_default.o: ../../lib/mod/heat_const.mod
heat_gill.o: ../../lib/mod/heat_const.mod
heat_gotm.o: ../../lib/mod/heat_const.mod
heat_temp.o: ../../lib/mod/heat_const.mod
hydrodynamic.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod ../../lib/mod/mkonst.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_fluidmud.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_nudging.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/mod_waves.mod \
		../../lib/mod/mod_zeta_system.mod \
		../../lib/mod/pkonst.mod ../../lib/mod/shympi.mod \
		../../lib/mod/shympi_debug.mod \
		../../lib/mod/tide.mod 
ice.o: ../../lib/mod/basin.mod
initialize.o: ../../lib/mod/basin.mod \
		../../lib/mod/coordinates.mod \
		../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod ../../lib/mod/mkonst.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/pkonst.mod ../../lib/mod/shympi.mod 
intp_tst.o: ../../lib/mod/basin.mod ../../lib/mod/mod_hydro.mod
iostr.o: ../../lib/mod/basin.mod ../../lib/mod/befor_after.mod \
		../../lib/mod/levels.mod ../../lib/mod/nls.mod \
		../../lib/mod/shympi.mod ../../lib/mod/simul.mod 
lagrange_back.o: ../../lib/mod/basin.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_connect.o: ../../lib/mod/basin.mod \
		../../lib/mod/femtime.mod \
		../../lib/mod/mod_debug.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_cont.o: ../../lib/mod/basin.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_decay.o: ../../lib/mod/femtime.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_dif.o: ../../lib/mod/basin.mod ../../lib/mod/femtime.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_flux.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_lagrange.mod \
		../../lib/mod/mod_layer_thickness.mod 
lagrange_init.o: ../../lib/mod/basin.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_inout.o: ../../lib/mod/basin.mod \
		../../lib/mod/femtime.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_lagrange.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/shyfile.mod 
lagrange_larve.o: ../../lib/mod/mod_lagrange.mod
lagrange_main.o: ../../lib/mod/basin.mod ../../lib/mod/femtime.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/lgr_sedim_module.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_lagrange.mod \
		../../lib/mod/shyfile.mod ../../lib/mod/shympi.mod 
lagrange_sedim.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_bstress.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_lagrange.mod \
		../../lib/mod/mod_sediment.mod 
lagrange_stk.o: ../../lib/mod/basin.mod \
		../../lib/mod/mod_lagrange.mod \
		../../lib/mod/mod_meteo.mod 
lagrange_track.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_lagrange.mod 
lagrange_util.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/mod_geom.mod 
lagrange_util_tr.o: ../../lib/mod/basin.mod
lagrange_utils.o: ../../lib/mod/mod_lagrange.mod
lagrange_vertical.o: ../../lib/mod/basin.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_lagrange.mod 
line_admin.o: ../../lib/mod/basin.mod ../../lib/mod/mod_geom.mod
lnkmpi.o: ../../lib/mod/basin.mod ../../lib/mod/mod_geom.mod \
		../../lib/mod/shympi.mod 
meteo_admin.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/meteo_forcing_module.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod 
meteo_forcing.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/shyfem_strings.mod \
		../../lib/mod/shympi.mod 
mkdistance.o: ../../lib/mod/basin.mod ../../lib/mod/mod_geom.mod \
		../../lib/mod/shympi.mod 
mod_petsc_global.o: ../../lib/mod/basin.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/petscdm.mod \
		../../lib/mod/petscdmlabel.mod \
		../../lib/mod/shympi.mod 
mod_petsc_system.o: ../../lib/mod/mod_petsc_global.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/petscdm.mod \
		../../lib/mod/petscdmlabel.mod \
		../../lib/mod/petscksp.mod \
		../../lib/mod/petscmat.mod \
		../../lib/mod/petscvec.mod \
		../../lib/mod/shympi.mod 
model3d_util.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/pkonst.mod \
		../../lib/mod/shympi.mod 
mpi_assert.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod 
mpi_tripple.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/shympi.mod 
mud_admin.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_fluidmud.mod 
nohydro.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/femtime.mod ../../lib/mod/levels.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_nohyd.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/mod_ts.mod \
		../../lib/mod/mod_zeta_system.mod \
		../../lib/mod/pkonst.mod 
nudge.o: ../../lib/mod/basin.mod ../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_nudging.mod \
		../../lib/mod/shympi.mod 
offline_data.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_offline.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod 
ous_admin.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/shyfile.mod ../../lib/mod/shympi.mod \
		../../lib/mod/simul.mod 
outputd.o: ../../lib/mod/femtime.mod
poisson.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod ../../lib/mod/mkonst.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/mod_zeta_system.mod \
		../../lib/mod/pkonst.mod 
residual.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/simul.mod 
restart.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/shympi.mod 
sedim_admin.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod ../../lib/mod/mkonst.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_bstress.mod \
		../../lib/mod/mod_debug.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/mod_sedtrans05.mod \
		../../lib/mod/mod_ts.mod \
		../../lib/mod/mod_waves.mod ../../lib/mod/para.mod 
shyfem.o: ../../lib/mod/basin.mod ../../lib/mod/befor_after.mod \
		../../lib/mod/clo.mod \
		../../lib/mod/coordinates.mod \
		../../lib/mod/evgeom.mod \
		../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_bclfix.mod \
		../../lib/mod/mod_bfm.mod \
		../../lib/mod/mod_bnd.mod \
		../../lib/mod/mod_bnd_aux.mod \
		../../lib/mod/mod_bndo.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_bound_geom.mod \
		../../lib/mod/mod_bstress.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_debug.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_aux.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_gotm_aux.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_keps.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_nohyd.mod \
		../../lib/mod/mod_nudging.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/mod_sediment.mod \
		../../lib/mod/mod_sinking.mod \
		../../lib/mod/mod_subset.mod \
		../../lib/mod/mod_test_zeta.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/mod_tvd.mod \
		../../lib/mod/mod_waves.mod \
		../../lib/mod/mod_zeta_system.mod \
		../../lib/mod/projection.mod \
		../../lib/mod/shympi.mod \
		../../lib/mod/shympi_debug.mod \
		../../lib/mod/subww3.mod ../../lib/mod/tide.mod 
simple_sedi.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_layer_thickness.mod 
simsys_lp.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_system.mod 
simsys_petsc.o: ../../lib/mod/mod_petsc_global.mod \
		../../lib/mod/mod_petsc_system.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/mod_system_interface.mod \
		../../lib/mod/petscdm.mod \
		../../lib/mod/petscdmlabel.mod \
		../../lib/mod/petscksp.mod \
		../../lib/mod/petscmat.mod \
		../../lib/mod/petscvec.mod \
		../../lib/mod/shympi.mod \
		../../lib/mod/shympi_time.mod 
simsys_spkaux.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/mod_system_interface.mod \
		../../lib/mod/shympi.mod 
sparskit_admin.o: ../../lib/mod/mod_system.mod \
		../../lib/mod/mod_system_interface.mod 
stability.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/pkonst.mod ../../lib/mod/shympi.mod 
subset.o: ../../lib/mod/basin.mod
substress.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_parwaves.mod \
		../../lib/mod/mod_ts.mod \
		../../lib/mod/mod_waves.mod \
		../../lib/mod/pkonst.mod ../../lib/mod/shympi.mod 
system_spk.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_system.mod \
		../../lib/mod/mod_system_interface.mod \
		../../lib/mod/shympi.mod 
test_zeta.o: ../../lib/mod/basin.mod ../../lib/mod/mod_hydro.mod \
		../../lib/mod/shympi.mod 
tidef.o: ../../lib/mod/basin.mod ../../lib/mod/coordinates.mod \
		../../lib/mod/iso8601.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_hydro.mod ../../lib/mod/tide.mod 
time_admin.o: ../../lib/mod/femtime.mod ../../lib/mod/shympi.mod
topological_admin.o: ../../lib/mod/basin.mod \
		../../lib/mod/evgeom.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/shympi.mod 
toxical.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_diff_visc_fric.mod 
trace.o: ../../lib/mod/basin.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_vel.mod 
transforms.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/shympi.mod 
tsfile_admin.o: ../../lib/mod/intp_fem_file.mod \
		../../lib/mod/levels.mod 
turbine.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_bound_dynamic.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_nudging.mod 
turbulence.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/femtime.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_area.mod \
		../../lib/mod/mod_diff_visc_fric.mod \
		../../lib/mod/mod_gotm_aux.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_keps.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/mod_ts.mod ../../lib/mod/pkonst.mod \
		../../lib/mod/shympi.mod 
tvd_admin.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/levels.mod \
		../../lib/mod/mod_hydro_print.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_tvd.mod ../../lib/mod/shympi.mod 
water_ren_time.o: ../../lib/mod/basin.mod \
		../../lib/mod/custom_dates.mod \
		../../lib/mod/evgeom.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_conz.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/shympi.mod ../../lib/mod/simul.mod 
waves_admin.o: ../../lib/mod/basin.mod \
		../../lib/mod/coordinates.mod \
		../../lib/mod/evgeom.mod ../../lib/mod/levels.mod \
		../../lib/mod/mod_depth.mod \
		../../lib/mod/mod_geom.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_vel.mod \
		../../lib/mod/mod_internal.mod \
		../../lib/mod/mod_layer_thickness.mod \
		../../lib/mod/mod_meteo.mod \
		../../lib/mod/mod_roughness.mod \
		../../lib/mod/mod_waves.mod \
		../../lib/mod/pkonst.mod 
wetdry.o: ../../lib/mod/basin.mod ../../lib/mod/evgeom.mod \
		../../lib/mod/mkonst.mod \
		../../lib/mod/mod_geom_dynamic.mod \
		../../lib/mod/mod_hydro.mod \
		../../lib/mod/mod_hydro_baro.mod \
		../../lib/mod/shympi.mod 

