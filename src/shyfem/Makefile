#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

FEMDIR = ../..
include $(FEMDIR)/Rules.make
FEMBIN = $(FEMDIR)/bin

#--------------------	EXES	---------------------

EXEC = shyfem

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

fem: $(EXEC) nothing

$(EXEC): $(MODULES) $(OBJ) $(LIBFS) shyfem.o
	 $(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

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
	-rm -f $(EXEC)
	-rm -f $(BINDIR)/$(EXEC)

cleanall: clean
	-rm -f $(MODDIR)/*.mod
	-rm -f $(LIBDIR)/libshy.a 

update:
	cp $(EXEC) /home/georg/work/shyfem_repo/shyfem-mpi/bin
	cp $(EXEC) /home/georg/work/shyfem_repo/shyfem-mpi/fem3d
	cp $(EXEC) /home/georg/work/shyfem_repo/shyfem/bin
	cp $(EXEC) /home/georg/work/shyfem_repo/shyfem/fem3d

#------------------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

