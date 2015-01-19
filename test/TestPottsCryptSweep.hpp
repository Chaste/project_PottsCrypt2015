#ifndef TESTPOTTSCRYPT_HPP_
#define TESTPOTTSCRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "TransitCellProliferativeType.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "WntConcentration.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "SloughingCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"

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
        unsigned start_sim = 0;
        unsigned num_sims = 1; //10
     	double mid_time = 1; //100
     	double end_time = 2; //200


     	//double temp[15] = {0.00001, 0.000031623, 0.0001, 0.00031623, 0.001, 0.0031623, 0.01, 0.031623, 0.1, 0.31623, 1.0, 3.1623, 10.0, 31.623, 100};
     	//unsigned max_temp_index = 15;
     	double temp[4] = {0.00001, 0.001, 0.1, 10.0};
     	unsigned max_temp_index = 4;

     	//double num_sweeps[12] = {1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500};
     	//unsigned max_num_sweeps_index = 10;
     	double dt[6] = {1, 1.0/3.0, 1.0/10.0, 1.0/30.0, 1.0/100.0, 1.0/300.0};
     	unsigned max_dt_index = 6;


     	// Stuff to output the number of cells in the crypt to a dat file.
        /** Results file cell velocities. */
        out_stream p_cell_number_file;
        OutputFileHandler output_file_handler("PottsCryptSweeps/", false);
        p_cell_number_file = output_file_handler.OpenOutputFile("cellnumbers.dat");

     	double number_of_cells_in_middle = 0.0;
     	double number_of_cells_at_end = 0.0;

		// Loop over Temp
		for (unsigned temp_index=0;  temp_index < max_temp_index; temp_index++)
		{
			std::cout << "\n Temp " << temp[temp_index] << ", " << std::flush;

			// Loop over Num Sweeps
			for (unsigned dt_index=0;  dt_index < max_dt_index; dt_index++)
			{
				std::cout << "\n\tDt " << dt[dt_index] << "... " << std::flush;

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
					simulator.SetSamplingTimestepMultiple((unsigned)(100.0/dt[dt_index]));
					simulator.SetOutputCellVelocities(true);

					//Create output directory
					std::stringstream out;
					out << "/Temp_"<< temp[temp_index] << "/Dt_" << dt[dt_index] << "/RunIndex_" << index;
					std::string output_directory = "Potts/CylindricalCrypt/Sweeps/" +  out.str();
					simulator.SetOutputDirectory(output_directory);

					// Create cell killer and pass in to simulation
					MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_length));
					simulator.AddCellKiller(p_killer);

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
					number_of_cells_in_middle += simulator.rGetCellPopulation().GetNumRealCells();

					// Run simulation to end
                    simulator.SetEndTime(end_time);
					simulator.Solve();
					number_of_cells_at_end += simulator.rGetCellPopulation().GetNumRealCells();

					// Tidy up
					WntConcentration<2>::Destroy();
					SimulationTime::Destroy();
					SimulationTime::Instance()->SetStartTime(0.0);
					RandomNumberGenerator::Destroy();
				}

			    number_of_cells_in_middle  /= (double) num_sims;
			    number_of_cells_at_end  /= (double) num_sims;

				*p_cell_number_file << temp[temp_index] << "\t" << dt[dt_index] << "\t" << number_of_cells_in_middle << "\t" << number_of_cells_at_end << "\n";

				number_of_cells_in_middle = 0.0;
			    number_of_cells_at_end = 0.0;
			}
        }
        std::cout << "\n" << std::flush;
        p_cell_number_file->close();
    }
};

#endif /*TESTPOTTSCRYPT_HPP_*/
