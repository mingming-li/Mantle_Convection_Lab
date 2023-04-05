void terminate();
void sync_processors();
void read_input_file();
void setup_mesh();
void setup_parser();

void get_time_independent_variables();
void initial_values();
void build_quad_tree();

struct TREE_NODE *new_tree_node();
void split_mul_lev();
void print_binary();
int morton();
void demorton();
void demorton_no_lev();
void printPostorder();
void get_Mcode_from_tree();
int get_Mcode_on_level();
void output_coord();
void construct_ien();
void construct_lm();
void construct_id();

void construct_v_lm();
void construct_v_id();
void construct_T_lm();
void construct_T_id();

void printCurrentLevel();
void boundary_condition();
void temperature_bc();
void velocity_bc();
void node_index_from_ien();
int append_level();
int append_leaf();
void sort();
void swap();
void swap_d();
void distribute_Mcode_among_processors();
void distribute_Mcode_evenly();
int get_Mcode_location();
int get_Mcode_level();
int get_Mcode_leaf();
int remove_dup();

void assign_Mcode_with_ip();
void move_Mcode_around();
void get_new_ip_nel();
void setup_processor_connection();

void global_index_Mcode();

int input_int();
int input_double();
int input_string();
int input_double_vector();

void construct_shape_function();
void shape_function();
void dshape_function_dxi_deta();
void dx_dxi_dz_deta();
void jacobi_in_integration();
double determinant();
void node_position();

void debug();
void initial_temperature();

void periodic_T();
void apply_T_bc_values();

void buoyancy_force();

void get_node_horiz_avg();
void remove_node_horiz_avg();

int cmp();

double Na();
int node_index();

void node_hashtable();

int hash_code();
struct HASH_ITEM *search();
void insert_hash_item();

void viscosity();

void velocity_pressure_KGF_matrix();
void K_matrix();
void F_matrix();
void G_matrix();
void H_matrix();

void get_element_K();
void get_element_F();
void get_element_G();
void get_element_H();

void solve_velocity_pressure();
void citcom_solver();
void initial_d_p();

void K_prod_d();
void G_prod_p();
void GT_prod_d();
void jacobi();
void multigrid();
void initial_error();
void MG();
int check_convergency();
void jacobi_smoother();
void interp_error();
void project_res();
void get_residual();
void link_eqs_across_levels();
void link_eq_ones();
void link_eq_twos();
void link_eq_fours();

void range_of_Mcode();
int what_ip_node_in();
void what_ip_el_share();
int max();
int min();
void link_eq_ones_local();
void link_eq_ones_remote();
void get_share_ip_nodes();
void link_proc_for_send();
void link_proc_for_recv();
void send_node_info();
void recv_node_info();
void remote_map_ones();
void mark_hangling_nodes();

void other_veq_contribution();
void other_eq_contribution();
void other_eq_contribution_old();
void other_eq_contribution_slow();
void other_lg_eq_contribution();

void find_sharing_nodes_eqs();

void assign_node_type();
void node_type();
void find_hangling_nodes();

void construct_connection();
void find_ip();
void proc_conn();
int locate_ip_in_nshareip();
void find_sharing_nodes();
void find_sharing_eqs();
void remote_KF_contribution();
int shared_eq();
void other_vel_contribution();
void what_elements_contain_eq();
void v_to_nodes();
void p_to_nodes();

void output_temperature();
void output();
void output_velocity();
int folder_exist();

void build_interp_conn();
void build_project_conn();

void release_memory();

void constrain_hangling_v_eq();
void hangling_Cv();
void build_Kf_matrix();
void enforce_constrain_eq();
void update_temperature();
void calculate_artificial_diffusivity();
void K_matrix_heat();
void F_matrix_heat();
void C_matrix_heat();
void get_element_K_heat();
void get_element_F_heat();
void get_element_C_heat();
void initial_T_dot();
void K_prod_d_heat();
void other_Teq_contribution();
void other_lg_Teq_contribution();
void enforce_constrain_Teq();
void jacobi_heat();

double advection_dt();
void Gauss_Seidel();
void Conj_Grad();
void Conj_Grad_heat();

void other_KFC_heat_contribution();
int shared_Teq();
void initial_data_structure();
void use_global_index_for_Mcode();
void project_Mcode();
void split_Mcode();
void initial_refine();
void append_level_to_Mcode();
int Mcode_lev();
int Mcode_loc();
int log4();

void adaptive_mesh_refine();

int find_element_ip();
int value_in_array();
int value_in_sorted_array();
void construct_ntype();
void elements_for_sharing_nodes();
void elements_for_sharing_eqs();
int find_elements_for_node();

void random_T_perturbation();
void periodic_T_with_random_perturb();
void gauss_seidel();
void conductive_T();
void setup_MPI_comm_vframe();
void setup_MPI_comm_Tframe();
void setup_project_comm_frame();
void setup_interp_comm_frame();
void assign_matrix_memory();
void interpolate_T_from_file();
void build_lagrange_K_matrix();
void add_lagrange_neq();
double global_v_prod();
double global_T_prod();
void find_nshare_for_eqs();
void find_nshare_for_Teqs();
void find_nshare_for_lg_eqs();
void find_nshare_for_lg_T_eqs();
void append_sharing_lg_eqs();
void build_lagrange_heat_KC_matrix();

double global_p_prod();

void visc_from_grain_size();
void compute_stress();
void compute_grain_size();

void element_area();

void what_elements_node_in();
void output_stress();
void output_viscosity();
void output_grain_size();
