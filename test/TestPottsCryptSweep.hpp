#ifndef TESTPOTTSCRYPT_HPP_
#define TESTPOTTSCRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "TransitCellProliferativeType.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "WntConcentration.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SloughingCellKiller.hpp"
#include "OnLatticeSimulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellIdWriter.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

class TestPottsCrypt : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:


    void TestMultiplePottsCrypt() throw (Exception)
    {
    	// Change to get different random seeds
        unsigned start_sim = 1;
        unsigned num_sims = 1; //10

     	double mid_time = 25;
     	double end_time = 50;


     	double cuberoot10 = 2.15443469;

     	double temp[10] = {0.001, 0.001*cuberoot10, 0.001*cuberoot10*cuberoot10, 0.01, 0.01*cuberoot10, 0.01*cuberoot10*cuberoot10, 0.1, 0.1*cuberoot10, 0.1*cuberoot10*cuberoot10, 1.0};// 3.1623, 10.0}; // 31.623, 100};
        unsigned max_temp_index = 10;

     	double dt[7] = {0.1, 0.01*cuberoot10*cuberoot10, 0.01*cuberoot10,  0.01, 0.001*cuberoot10*cuberoot10, 0.001*cuberoot10, 0.001};
        unsigned max_dt_index = 7;

     	// Code to output the number of cells in the crypt to a dat file.
        /** Results file cell velocities. */
        out_stream p_cell_number_file;
        OutputFileHandler output_file_handler("Potts/CylindricalCrypt/Sweeps/", false);
        p_cell_number_file = output_file_handler.OpenOutputFile("cellnumbers.dat");

     	// Loop over dt
		for (unsigned dt_index=0;  dt_index < max_dt_index;  dt_index++)
		{
			std::cout << "\nDt " << dt[dt_index] << "... " << std::flush;

			// Loop over Temp
			for (unsigned temp_index=0;  temp_index < max_temp_index; temp_index++)
			{
				std::cout << "\n\tTemp " << temp[temp_index] << ", " << std::flush;

		     	double number_of_cells_in_middle = 0.0;
		     	double number_of_cells_at_end = 0.0;

		     	// Loop over seed
			    for(unsigned index=start_sim; index < start_sim + num_sims; index++)
				{
					std::cout << " Run number " << index << "... " << std::flush;


					// Re seed the random number generator
					RandomNumberGenerator::Instance()->Reseed(100*index);

					double crypt_length = 100;

					// Create a simple 2D PottsMesh
					PottsMeshGenerator<2> generator(50, 10, 5, 110, 20, 5, 1, 1, 1, true, true);
					PottsMesh<2>* p_mesh = generator.GetMesh();

                    MAKE_PTR(TransitCellProliferativeType, p_transit_type);

                    // Create cells
					std::vector<CellPtr> cells;
					CellsGenerator<SimpleWntCellCycleModel, 2> cells_generator;
					cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);


					// Alter cells properties
					for (unsigned i=0; i<cells.size(); i++)
					{
						dynamic_cast<SimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetTransitCellG1Duration(6);
						dynamic_cast<SimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetWntTransitThreshold(2.0/3.0);
					}

					// Create cell population
					PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
					cell_population.SetNumSweepsPerTimestep(1);
					cell_population.SetTemperature(temp[temp_index]);

			        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
			        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
			        cell_population.AddCellWriter<CellMutationStatesWriter>();
			        cell_population.AddCellWriter<CellVolumesWriter>();
					cell_population.AddCellWriter<CellIdWriter>();

					// Create an instance of a Wnt concentration
					WntConcentration<2>::Instance()->SetType(LINEAR);
					WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
					WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

					// Set up cell-based simulation
					OnLatticeSimulation<2> simulator(cell_population);
					simulator.SetDt(dt[dt_index]);
					simulator.SetSamplingTimestepMultiple((unsigned)(1.0/dt[dt_index]));
					simulator.SetOutputCellVelocities(true);

					//Create output directory
					std::stringstream out;
					out <<  "/Dt_" << dt[dt_index] << "/Temp_"<< temp[temp_index] << "/RunIndex_" << index;
					std::string output_directory = "Potts/CylindricalCrypt/Sweeps/" +  out.str();
					simulator.SetOutputDirectory(output_directory);

					// Create cell killer and pass in to simulation
					MAKE_PTR_ARGS(SloughingCellKiller<2>, p_sloughing_killer, (&cell_population, crypt_length));
					simulator.AddCellKiller(p_sloughing_killer);

					// Create update rules and pass to the simulation
					MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
					p_volume_constraint_update_rule->SetMatureCellTargetVolume(25);
					p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.1); //Default is 0.5
					simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
					MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
					simulator.AddPottsUpdateRule(p_adhesion_update_rule);

					// Run simulation to middle
					simulator.SetEndTime(mid_time);
					simulator.Solve();

					// Get number of elements of non zero size
					PottsMesh<2>& potts_mesh = static_cast<PottsMesh<2>&>(simulator.rGetCellPopulation().rGetMesh());
					unsigned local_num_cells_in_middle=0;
					for (PottsMesh<2>::PottsElementIterator elem_iter = potts_mesh.GetElementIteratorBegin();
					   elem_iter != potts_mesh.GetElementIteratorEnd();
					   ++elem_iter)
					{
						if (elem_iter->GetNumNodes()>0)
						{
							local_num_cells_in_middle ++;
						}
					}
					number_of_cells_in_middle += local_num_cells_in_middle;

					// Run simulation to end
                    simulator.SetEndTime(end_time);
					simulator.Solve();


					// Get number of elements of non zero size
					PottsMesh<2>& potts_mesh_end = static_cast<PottsMesh<2>&>(simulator.rGetCellPopulation().rGetMesh());

					unsigned local_num_cells_at_end=0;
					for (PottsMesh<2>::PottsElementIterator elem_iter = potts_mesh_end.GetElementIteratorBegin();
					   elem_iter != potts_mesh_end.GetElementIteratorEnd();
					   ++elem_iter)
					{
						if (elem_iter->GetNumNodes()>0)
						{
							local_num_cells_at_end ++;
						}
					}
					number_of_cells_at_end += local_num_cells_at_end;

					// Tidy up
					WntConcentration<2>::Destroy();
					SimulationTime::Destroy();
					SimulationTime::Instance()->SetStartTime(0.0);
					RandomNumberGenerator::Destroy();

					*p_cell_number_file << temp[temp_index] << "\t" << dt[dt_index] << "\t" << index << "\t" << local_num_cells_in_middle << "\t" << local_num_cells_at_end << "\n";

				}
				number_of_cells_in_middle = 0.0;
			    number_of_cells_at_end = 0.0;
			}
        }
        std::cout << "\n" << std::flush;
        p_cell_number_file->close();
    }
};

#endif /*TESTPOTTSCRYPT_HPP_*/
