/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	//cell_defaults.functions.update_phenotype = phenotype_function; 
	//cell_defaults.functions.custom_cell_rule = custom_function; 
	//cell_defaults.functions.contact_function = contact_function; 
	
	Cell_Definition* pCD = find_cell_definition("epithelial");

	pCD->functions.update_phenotype = NULL;
	pCD->functions.custom_cell_rule = custom_function;
	pCD->functions.contact_function = contact_function;
	
	pCD->phenotype.mechanics.attachment_elastic_constant = 0.03;

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	/*Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return;*/

	Cell* pCell;
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.8 * 2.0 * cell_radius; 
	double initial_tumor_radius = 100;

	std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius,initial_tumor_radius);
    
    std::cout << "NUMBER OF CELLS : " << positions.size() << " __________" << std::endl;
    for( int i=0; i < positions.size(); i++ )
    {
        pCell = create_cell(get_cell_definition("epithelial"));
		int ext_delta = pCell->custom_data.find_variable_index( "ext_delta" );
        int ext_notch = pCell->custom_data.find_variable_index( "ext_notch" );
        int ext_jagged = pCell->custom_data.find_variable_index( "ext_jagged" );
		int delta_avail_attach = pCell->custom_data.find_variable_index( "delta_avail_attach" );
        int notch_avail_attach = pCell->custom_data.find_variable_index( "notch_avail_attach" );
    	int jagged_avail_attach = pCell->custom_data.find_variable_index( "jagged_avail_attach" );
        int delta = pCell->custom_data.find_variable_index( "delta" );
        int notch = pCell->custom_data.find_variable_index( "notch" );
        int jagged = pCell->custom_data.find_variable_index( "jagged" );
        int NICD = pCell->custom_data.find_variable_index( "NICD" );
        int mig_speed = pCell->custom_data.find_variable_index( "mig_speed" ); 
		int pers_time = pCell->custom_data.find_variable_index( "pers_time" );
        pCell->assign_position( positions[i] );

        pCell->custom_data[delta] = (int) (UniformRandom()*1500);
        pCell->custom_data[notch] = (int) (UniformRandom()*1500);
        pCell->custom_data[jagged] = (int) (UniformRandom()*1500);
        pCell->custom_data[NICD] = (int) (UniformRandom()*200);
		pCell->custom_data[mig_speed] = 1;
		pCell->custom_data[pers_time] = 0.5;

		//Inititalizing _avail_attach and ext_
			//no cell attached at t = 0
		pCell->custom_data[delta_avail_attach] = 0;	
		pCell->custom_data[notch_avail_attach] = 0;
		pCell->custom_data[jagged_avail_attach] = 0;
		pCell->custom_data[ext_delta] = 0;
        pCell->custom_data[ext_notch] = 0;
        pCell->custom_data[ext_jagged] = 0;
        
        /*double cell_volume = pCell->phenotype.volume.total;
        pCell->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index]= pCell->custom_data[i_Oxy_i] * cell_volume;
        pCell->phenotype.molecular.internalized_total_substrates[glucose_substrate_index]= pCell->custom_data[i_Glu_i] * cell_volume;
        pCell->phenotype.molecular.internalized_total_substrates[lactate_substrate_index]= pCell->custom_data[i_Lac_i] * cell_volume;*/
		
        pCell->phenotype.intracellular->start();
		pCell->phenotype.intracellular->set_parameter_value("D",pCell->custom_data[delta]);
		pCell->phenotype.intracellular->set_parameter_value("N",pCell->custom_data[notch]);
		pCell->phenotype.intracellular->set_parameter_value("J",pCell->custom_data[jagged]);
        pCell->phenotype.intracellular->set_parameter_value("I",pCell->custom_data[NICD]);	
		pCell->phenotype.intracellular->set_parameter_value("Dext",pCell->custom_data[ext_delta]);
		pCell->phenotype.intracellular->set_parameter_value("Next",pCell->custom_data[ext_notch]);
		pCell->phenotype.intracellular->set_parameter_value("Jext",pCell->custom_data[ext_jagged]);
		pCell->phenotype.intracellular->set_parameter_value("mig_speed",pCell->custom_data[mig_speed]); 
		pCell->phenotype.intracellular->set_parameter_value("pers_time",pCell->custom_data[pers_time]);
       
    }
    
	return;  
}


void update_intracellular()
{	
    // BioFVM Indices
    //static int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );
    //static int glucose_substrate_index = microenvironment.find_density_index( "glucose" ); 
    //static int lactate_substrate_index = microenvironment.find_density_index( "lactate");

    #pragma omp parallel for   
    for( int i=0; i < (*all_cells).size(); i++ )
    {	
		double ext_delta_val = 0;
		double ext_notch_val = 0;
		double ext_jagged_val = 0;

		//Get the indices for custom cell data
        static int ext_delta = (*all_cells)[i]->custom_data.find_variable_index( "ext_delta" );
        static int ext_notch = (*all_cells)[i]->custom_data.find_variable_index( "ext_notch" );
        static int ext_jagged = (*all_cells)[i]->custom_data.find_variable_index( "ext_jagged" );
		static int delta_avail_attach_me = (*all_cells)[i]->custom_data.find_variable_index( "delta_avail_attach" );
        static int notch_avail_attach_me = (*all_cells)[i]->custom_data.find_variable_index( "notch_avail_attach" );
        static int jagged_avail_attach_me = (*all_cells)[i]->custom_data.find_variable_index( "jagged_avail_attach" );
        static int delta = (*all_cells)[i]->custom_data.find_variable_index( "delta" );
        static int notch = (*all_cells)[i]->custom_data.find_variable_index( "notch" );
        static int jagged = (*all_cells)[i]->custom_data.find_variable_index( "jagged" );
		static int NICD = (*all_cells)[i]->custom_data.find_variable_index( "NICD" );
		static int mig_speed = (*all_cells)[i]->custom_data.find_variable_index( "mig_speed" ); 
		static int pers_time = (*all_cells)[i]->custom_data.find_variable_index( "pers_time" );
		
		//Get the list of cell attached 
		std::vector<Cell*> cell_attached = (*all_cells)[i]->state.attached_cells;
		//std::cout<<" Check 5"<<std::endl;

		//Calculate the ext_notch, ext_delta and ext_jagged
		for( int j=0; j < cell_attached.size(); j++)
		{	
			int delta_avail_attach_other = cell_attached[j]->custom_data.find_variable_index( "delta_avail_attach" );
        	int notch_avail_attach_other = cell_attached[j]->custom_data.find_variable_index( "notch_avail_attach" );
        	int jagged_avail_attach_other = cell_attached[j]->custom_data.find_variable_index( "jagged_avail_attach" );
			ext_delta_val += cell_attached[j]->custom_data[delta_avail_attach_other];
			ext_notch_val += cell_attached[j]->custom_data[notch_avail_attach_other];
			ext_jagged_val += cell_attached[j]->custom_data[jagged_avail_attach_other];	
		}

        //std::cout<<" Check 6"<<std::endl;
        if( (*all_cells)[i]->is_out_of_domain == false  )
        {
            // Cell Volume
            //double cell_volume = (*all_cells)[i]->phenotype.volume.total;
            
            // Intracellular Concentrations
            double delta_val_int = (*all_cells)[i]->custom_data[delta];
            double notch_val_int = (*all_cells)[i]->custom_data[notch];
            double jagged_val_int = (*all_cells)[i]->custom_data[jagged];
			double NICD_val_int = (*all_cells)[i]->custom_data[NICD];
            
            // Update SBML 
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("Dext",ext_delta_val);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Next",ext_notch_val);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Jext",ext_jagged_val);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("D",delta_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("N",notch_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("J",jagged_val_int);
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("mig_speed",(*all_cells)[i]->custom_data[mig_speed]); 
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("pers_time",(*all_cells)[i]->custom_data[pers_time]);
                        
            // SBML Simulation
            (*all_cells)[i]->phenotype.intracellular->update();
            // Phenotype Simulation
            (*all_cells)[i]->phenotype.intracellular->update_phenotype_parameters((*all_cells)[i]->phenotype);
            
            //std::cout << "Before Intracellular Oxygen : " <<(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index]/cell_volume << std::endl;
            
            // Internalized Chemical Update After SBML Simulation
            //(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen") * cell_volume;
            //(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glucose_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose") * cell_volume;
            //(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lactate_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate") * cell_volume;
            
            //std::cout << "SBML Energy : " <<(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") << std::endl;
            /* if ( (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") >100 )
            {
                std::cout << "SBML Energy : " <<(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") << "  - Cell position : " << (*all_cells)[i]->position << std::endl;
            } */
            
            //Save custom data
            (*all_cells)[i]->custom_data[delta] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("D");
            (*all_cells)[i]->custom_data[notch] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("N");
            (*all_cells)[i]->custom_data[jagged] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("J");
            (*all_cells)[i]->custom_data[NICD] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("I");

			//updating the ext_delta custom cell data parameters
			(*all_cells)[i]->custom_data[ext_delta] = ext_delta_val;
			(*all_cells)[i]->custom_data[ext_notch] = ext_notch_val;
			(*all_cells)[i]->custom_data[ext_jagged] = ext_jagged_val;

			//updating _avail_attach receptors 
			//(*all_cells)[i]->custom_data[delta_avail_attach_me] = (*all_cells)[i]->custom_data[delta]/((*all_cells)[i]->state.number_of_attached_cells()+1);
			//(*all_cells)[i]->custom_data[notch_avail_attach_me] = (*all_cells)[i]->custom_data[notch]/((*all_cells)[i]->state.number_of_attached_cells()+1);
			//(*all_cells)[i]->custom_data[jagged_avail_attach_me] = (*all_cells)[i]->custom_data[jagged]/((*all_cells)[i]->state.number_of_attached_cells()+1);
			update_avail_attach((*all_cells)[i]);
		}
    }

	return;
}

std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
			tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=0;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
		}
	}
	return cells;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ //return paint_by_number_cell_coloring(pCell); 
	
	//get the index for notch cell data
	static int nicd = pCell->custom_data.find_variable_index( "NICD" );
	static int low_nicd = pCell->custom_data.find_variable_index( "low_nicd" );
	static int high_nicd = pCell->custom_data.find_variable_index( "high_nicd" );

	//color based in the value of notch
	if(pCell->custom_data[nicd] < pCell->custom_data[low_nicd])
	{	return{"blue","black","blue","blue"};	}

	if(pCell->custom_data[nicd] > pCell->custom_data[low_nicd] && pCell->custom_data[nicd] < pCell->custom_data[high_nicd])
	{	return{"Red","black","Red","Red"};	}

	if(pCell->custom_data[nicd] > pCell->custom_data[high_nicd])
	{	return{"yellow","black","yellow","yellow"};	}	

	return{"grey","black","grey","grey"};
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 	
	//Attachment
	
	//list of nearby cells
	std::vector<Cell*> nearby = pCell->nearby_interacting_cells();

	//value of nicd for the cell
	static int nicd = pCell->custom_data.find_variable_index("NICD");
	static int low_nicd = pCell->custom_data.find_variable_index( "low_nicd" );
	static int high_nicd = pCell->custom_data.find_variable_index( "high_nicd" );

	double nicd_val = pCell->custom_data[nicd];
	int n = 0;
	double check = 0.0;
	//std::cout<<"Near by "<<nearby.size()<<std::endl;
	
	while (n < nearby.size())
	{	
		if(nearby[n] != pCell)// && pCell->state.number_of_attached_cells() < 7) 
		{	
			//std::cout<<"Check 1"<<std::endl;
			//pCell->attach_cell(nearby[n]);
		
			//value of nicd for the nearby cell
			int nicd_neigh = nearby[n]->custom_data.find_variable_index("NICD");
			double nicd_neigh_val = nearby[n]->custom_data[nicd_neigh]; 

			//pCell -> epithelial
			if( nicd_val <= pCell->custom_data[low_nicd])		
			{	
				//nearby is epithelial
				if(nicd_neigh_val <= nearby[n]->custom_data[low_nicd])		
				{	
					pCell->attach_cell(nearby[n]);
					update_avail_attach(pCell);
					update_avail_attach(nearby[n]);
				}
				//nearby is hybrid
				else if(nicd_neigh_val > nearby[n]->custom_data[low_nicd] && nicd_neigh_val <= nearby[n]->custom_data[high_nicd])
				{	check = UniformRandom();	
					if(check <= 0.75)	
					{	
						pCell->attach_cell(nearby[n]);	
						nearby[n]->phenotype.motility.migration_speed = 0;
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
				
				//nearby is mesenchymal
				else if(nicd_neigh_val > nearby[n]->custom_data[high_nicd])
				{	check = UniformRandom();	
					if(check <= 0.6)
					{	
						pCell->attach_cell(nearby[n]);
						nearby[n]->phenotype.motility.migration_speed = 0;
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
			}

			//pCell -> hybrid
			else if(nicd_val > pCell->custom_data[low_nicd] && nicd_val <= pCell->custom_data[high_nicd])
			{
				if(nicd_neigh_val <= nearby[n]->custom_data[low_nicd])		
				{	check = UniformRandom();	
					if(check <= 0.75)
					{	
						pCell->attach_cell(nearby[n]);
						pCell->phenotype.motility.migration_speed = 0;
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
				else if(nicd_neigh_val > nearby[n]->custom_data[low_nicd] && nicd_neigh_val <= nearby[n]->custom_data[high_nicd])
				{	check = UniformRandom();	
					if(check <= 0.75)
					{	
						pCell->attach_cell(nearby[n]);	
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
				else if(nicd_neigh_val > nearby[n]->custom_data[high_nicd])
				{	check = UniformRandom();	
					if(check <= 0.6)
					{	
						pCell->attach_cell(nearby[n]);
						pCell->phenotype.motility.migration_speed = 0.5;
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
			}

			//pCell -> mesenchymal
			else if( nicd_val > pCell->custom_data[high_nicd])
			{
				if(nicd_neigh_val<=nearby[n]->custom_data[low_nicd])		
				{	check = UniformRandom();	
					if(check <= 0.6)
					{	
						pCell->attach_cell(nearby[n]);
						pCell->phenotype.motility.migration_speed = 0;
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
				else if(nicd_neigh_val > nearby[n]->custom_data[low_nicd] && nicd_neigh_val <= nearby[n]->custom_data[high_nicd])
				{	check = UniformRandom();	
					if(check <= 0.6)
					{	
						pCell->attach_cell(nearby[n]);
						nearby[n]->phenotype.motility.migration_speed = 0.5;
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
				else if(nicd_neigh_val>nearby[n]->custom_data[low_nicd])
				{	check = UniformRandom();	
					if(check <= 0.3)
					{	
						pCell->attach_cell(nearby[n]);	
						update_avail_attach(pCell);
						update_avail_attach(nearby[n]);
					}
				}
			}
		}
		n++;
	}

	//Detachment

	//Get the list of cell attached 
	std::vector<Cell*> cell_attached = pCell->state.attached_cells;

	n = 0;
	check = 0.0;

	//std::cout<<"Cells attached "<<cell_attached.size()<<std::endl;

	while (n < cell_attached.size())
	{	//std::cout<<"Check 2"<<std::endl;
		if(cell_attached[n] != pCell)
		{	
			//value of nicd for the nearby cell
			int nicd_attach = cell_attached[n]->custom_data.find_variable_index("NICD");
			double nicd_attach_val = cell_attached[n]->custom_data[nicd_attach]; 

			//pCell -> epithelial
			if( nicd_val<=pCell->custom_data[low_nicd])		
			{
				if(nicd_attach_val<=cell_attached[n]->custom_data[low_nicd])		
				{	
					check = UniformRandom();	
					if(check <= 0.1)
					{	
						pCell->detach_cell(cell_attached[n]);
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
				else if(nicd_attach_val>cell_attached[n]->custom_data[low_nicd] && nicd_attach_val<=cell_attached[n]->custom_data[high_nicd])
				{	
					check = UniformRandom();	
					if(check <= 0.3)
					{	
						pCell->detach_cell(cell_attached[n]);	
						cell_attached[n]->phenotype.motility.migration_speed = 1;
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
				else if(nicd_attach_val>cell_attached[n]->custom_data[high_nicd])
				{	check = UniformRandom();	
					if(check <= 0.5)
					{	
						pCell->detach_cell(cell_attached[n]);	
						cell_attached[n]->phenotype.motility.migration_speed = 0.5;
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
			}

			//pCell -> hybrid
			else if( nicd_val > pCell->custom_data[low_nicd] && nicd_val <= pCell->custom_data[high_nicd])
			{
				if(nicd_attach_val<=cell_attached[n]->custom_data[low_nicd])		
				{	
					check = UniformRandom();	
					if(check <= 0.3)
					{	
						pCell->detach_cell(cell_attached[n]);
						pCell->phenotype.motility.migration_speed = 1;
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
				else if(nicd_attach_val>cell_attached[n]->custom_data[low_nicd] && nicd_attach_val<=cell_attached[n]->custom_data[high_nicd])
				{	
					check = UniformRandom();	
					if(check <= 0.3)
					{	
						pCell->detach_cell(cell_attached[n]);	
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
				else if(nicd_attach_val>cell_attached[n]->custom_data[high_nicd])
				{	
					check = UniformRandom();	
					if(check <= 0.5)
					{	
						pCell->detach_cell(cell_attached[n]);	
						pCell->phenotype.motility.migration_speed = 1;
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
			}

			//pCell -> mesenchymal
			else if( nicd_val > pCell->custom_data[high_nicd])
			{
				if(nicd_attach_val <= cell_attached[n]->custom_data[low_nicd])		
				{	
					check = UniformRandom();	
					if(check <= 0.5)
					{	
						pCell->detach_cell(cell_attached[n]);
						pCell->phenotype.motility.migration_speed = 0.5;
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
				else if(nicd_attach_val > cell_attached[n]->custom_data[low_nicd] && nicd_attach_val<=cell_attached[n]->custom_data[high_nicd])
				{	
					check = UniformRandom();	
					if(check <= 0.5)
					{	
						pCell->detach_cell(cell_attached[n]);	
						cell_attached[n]->phenotype.motility.migration_speed = 1;
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
				else if(nicd_attach_val>cell_attached[n]->custom_data[high_nicd])
				{	
					check = UniformRandom();	
					if(check <= 0.7)
					{	
						pCell->detach_cell(cell_attached[n]);	
						update_avail_attach(pCell);
						update_avail_attach(cell_attached[n]);
					}
				}
			}
		}
		n++;
	}

	return; 
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return standard_elastic_contact_function(pMe,phenoMe,pOther,phenoOther,dt); } 

//To update the _avail_attach values
void update_avail_attach(Cell* pCell)
{
	static int delta_avail_attach = pCell->custom_data.find_variable_index( "delta_avail_attach" );
	static int notch_avail_attach = pCell->custom_data.find_variable_index( "notch_avail_attach" );
	static int jagged_avail_attach = pCell->custom_data.find_variable_index( "jagged_avail_attach" );
	static int delta = pCell->custom_data.find_variable_index( "delta" );
	static int notch = pCell->custom_data.find_variable_index( "notch" );
	static int jagged = pCell->custom_data.find_variable_index( "jagged" );

	//updating _avail_attach receptors 
	if( pCell->state.number_of_attached_cells() > 0)
	{
		pCell->custom_data[delta_avail_attach] = pCell->custom_data[delta]/(pCell->state.number_of_attached_cells());
		pCell->custom_data[notch_avail_attach] = pCell->custom_data[notch]/(pCell->state.number_of_attached_cells());
		pCell->custom_data[jagged_avail_attach] = pCell->custom_data[jagged]/(pCell->state.number_of_attached_cells());
	}
	else
	{
		pCell->custom_data[delta_avail_attach] = 0;
		pCell->custom_data[notch_avail_attach] = 0;
		pCell->custom_data[jagged_avail_attach] = 0;
	}

}
