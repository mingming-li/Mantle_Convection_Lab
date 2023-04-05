#include "global_variables.h"

void read_input_file()
{
	char temp[1];

	input_string("data_dir",G.data_dir,".");
	input_string("dimension",G.dimension,"2d");
	input_string("geometry",G.geometry,"cart");
	input_string("solver",G.solver,"jc");

	input_int("level",&G.max_level,"1");
	input_int("mg_min_level",&G.min_level,"1");
	if(G.min_level >= G.max_level) 
	{
		fprintf(stderr,"min_level must be smaller than max_level\n");
		terminate();
	}

	input_int("nprox",&G.npx,"2");
	input_int("nproz",&G.npz,"2");

	input_double("model_length",&G.global.length,"2.0");
	input_double("model_height",&G.global.height,"1.0");

	input_string("top_v_bc",G.bc.top_v_bc,"sv");
	input_string("bot_v_bc",G.bc.bot_v_bc,"sv");
	input_string("lef_v_bc",G.bc.lef_v_bc,"vs");
	input_string("rig_v_bc",G.bc.rig_v_bc,"vs");
	input_double_vector("top_v_bc_val",NSD,G.bc.top_v_bc_val,"0.0,0.0");
	input_double_vector("bot_v_bc_val",NSD,G.bc.bot_v_bc_val,"0.0,0.0");
	input_double_vector("lef_v_bc_val",NSD,G.bc.lef_v_bc_val,"0.0,0.0");
	input_double_vector("rig_v_bc_val",NSD,G.bc.rig_v_bc_val,"0.0,0.0");

	input_string("top_T_bc",temp,"T"); G.bc.top_T_bc = temp[0];
	input_string("bot_T_bc",temp,"T"); G.bc.bot_T_bc = temp[0];
	input_string("lef_T_bc",temp,"F"); G.bc.lef_T_bc = temp[0];
	input_string("rig_T_bc",temp,"F"); G.bc.rig_T_bc = temp[0];

	input_double("top_T_bc_val",&G.bc.top_T_bc_val,"0.0");
	input_double("bot_T_bc_val",&G.bc.bot_T_bc_val,"1.0");
	input_double("lef_T_bc_val",&G.bc.lef_T_bc_val,"0.0");
	input_double("rig_T_bc_val",&G.bc.rig_T_bc_val,"0.0");

	input_int("T_ini",&G.control.T_ini,"1");
        input_double("T_mantle", &G.control.T_mantle,"0.5");
        input_double("T_mantle_perturb", &G.control.T_mantle_perturb,"0.01");
        input_int("T_wave_number", &G.control.T_wave_number,"6");


        input_double("Rayleigh", &G.Ra,"1e7");

        input_int("use_hash", &G.control.use_hash,"0");

	input_double("A", &G.viscosity.A, "0.0");
	input_double("refT", &G.viscosity.refT, "0.5");
	input_double("lower_mantle_jump", &G.viscosity.lm_jump, "1.0");
	input_double("highest_vis", &G.viscosity.highest, "1.0");
	input_double("lowest_vis", &G.viscosity.lowest, "1.0");

	input_int("use_dense_matrix",&G.control.use_dense_matrix,"1");


	input_int("mg_coarse_iteration",&G.control.mg_coarse_iteration,"4");
	input_int("mg_finest_iteration",&G.control.mg_finest_iteration,"4");

        input_int("max_timestep", &G.max_timestep,"0");
        input_int("save_spacing", &G.save_spacing,"1");

	input_double("Q",&G.control.Q,"0.0");
	input_double("theta",&G.control.theta,"0.0");

	input_int("AMR",&G.control.AMR,"0");

	input_int("mpi_KF_method",&G.control.mpi_KF_method,"2");

	input_int("cg_precondition",&G.control.cg_precondition,"0");
	input_int("Tcg_precondition",&G.control.Tcg_precondition,"0");

	input_double("accuracy",&G.control.accuracy,"0.001");
	input_int("auto_accuracy",&G.control.auto_accuracy,"0");
	input_double("citcom_solver_accu",&G.control.citcom_accu,"0.001");
	input_int("develop_mode",&G.control.dev_mode,"1");
	input_int("check_v_only",&G.control.check_v_only,"0");
	input_int("max_citcom_cycle",&G.control.max_citcom_cycle,"375");
	input_int("mg_cycle",&G.control.mg_cycle,"2");
	input_int("max_cg_cycle",&G.control.max_cg_cycle,"200");

	input_string("old_T_file",G.old_T_file,"");
	input_int("old_T_level",&G.old_T_level,"1");
	input_int("old_T_timestep",&G.old_T_timestep,"0");


	input_int("lg_mul",&G.lg_mul,"0");
	input_int("grain_damage",&G.control.grain_damage,"0");
	input_int("show_convg",&G.control.show_convg,"0");
	input_int("zero_P",&G.control.zero_P,"1");
	return;
}
