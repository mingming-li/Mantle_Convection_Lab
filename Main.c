#include "global_variables.h"

int main(int argc,char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &G.ip);
	MPI_Comm_size(MPI_COMM_WORLD, &G.np);

	setup_parser(argc,argv);
	read_input_file();

	if(G.np!=G.npx * G.npz)
        {
                fprintf(stderr,"Number of processors do not match the command line\n");
                terminate();
        }

	if(G.ip==0)
		start_time=MPI_Wtime();


	initial_data_structure();
	get_time_independent_variables();


	G.timestep = 0;
	G.solution_cycle = 0;
	while(G.timestep <= G.max_timestep)
	{
		if(G.ip == 0)
			fprintf(stderr,"\ntimestep = %d, CPU time = %.0fs, model time = %e\n",G.timestep, MPI_Wtime() - start_time, G.model_time);

		if(G.solution_cycle == 0 || G.control.AMR == 1)
		{
			build_quad_tree();
			if(G.ip == 0)
				fprintf(stderr, "Tree setup done!\n");
			distribute_Mcode_among_processors();
			construct_ien();
			node_index_from_ien();
			node_position();
			element_area();
			what_elements_node_in();
			boundary_condition();
			construct_id();
			construct_lm();
			if(G.ip == 0)
				fprintf(stderr, "Setup connection\n");
			construct_connection();
			construct_shape_function();
			assign_matrix_memory();
		}


		if(G.solution_cycle == 0)
			initial_temperature();

		buoyancy_force();
		viscosity();
		velocity_pressure_KGF_matrix();
		solve_velocity_pressure();
		if(G.timestep % G.save_spacing == 0)
			output();
		update_temperature();
		if(G.control.AMR == 1)
			release_memory();
		G.timestep++;
		G.solution_cycle++;
	}

	if(G.ip==0)
	{
		fprintf(stderr,"%f seconds of CPU time used in total\n",MPI_Wtime()-start_time);
		//printf("%d x %d\n",G.gnox,G.gnoz);
	}

	MPI_Finalize();
	return 0;
}

