#ifndef TESTPOTTSCRYPTMUTANT_HPP_
#define TESTPOTTSCRYPTMUTANT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "MutantBaseTrackerModifier.hpp"
#include "CellShapeOutputModifier.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "WntConcentration.hpp"
#include "FixedSimpleWntCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "MutantCellPottsUpdateRule.hpp"
#include "SloughingCellKiller.hpp"
#include "OnLatticeSimulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellVolumesWriter.hpp"

#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CommandLineArguments.hpp"
#include "Debug.hpp"

class TestPottsCryptMutant : public AbstractCellBasedTestSuite
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

    void TestMultipleMutantPottsCrypts() throw (Exception)
    {

    	// Uncoment these lines to make an executable with arguments, useful for sweeping

//        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-run_index"));
//        unsigned start_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-run_index");
//
//        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_runs"));
//        unsigned num_runs = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-num_runs");

    	unsigned start_index = 0;
    	unsigned num_runs = 1;

        std::string main_directory = "PottsCryptMutant/";
        std::string steady_state_output_directory, output_directory;

        unsigned crypt_length = 100; //100
        unsigned crypt_width = 50; //50
        unsigned element_size = 5;

        double time_to_steady_state = 50.0; //50
        double time_after_mutations = 50.0; //50

        double blob_radius = 10.0; // Lattice sites so 2 CDs

        unsigned num_blob_heights = 1;
        double blob_heights[2] = {20.0,60.0}; // Lattice sites so 4 and 12 CDs


        unsigned num_drag_ratios = 1;
        double drag_ratios[19] = {1.0, 1.5, 2.0, 2.5, 3.0,3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};

        for(unsigned index=start_index; index < start_index + num_runs; index++)
        {
            std::cout << "\nExperiment number " << index << "... " << std::flush;

            // loop over drag
            for (unsigned drag_index= 0; drag_index < num_drag_ratios; drag_index++)
            {
                std::cout << "Drag " << drag_ratios[drag_index] << ", " << std::flush;

                // loop over mutation blob heights
                for (unsigned height_index=0; height_index < num_blob_heights; height_index++)
                {
                    std::cout << "Height " << blob_heights[height_index] << "... " << std::flush;

                    // To be extra careful, we reseed the random number generator
                    RandomNumberGenerator::Instance()->Reseed(100*index);

                    // Create a simple 2D PottsMesh
                    PottsMeshGenerator<2> generator(crypt_width, crypt_width/element_size, element_size, crypt_length +10 , crypt_length/element_size, element_size, 1, 1, 1, true, true);
                    PottsMesh<2>* p_mesh = generator.GetMesh();

                    MAKE_PTR(TransitCellProliferativeType, p_transit_type);

                    // Create cells
                    std::vector<CellPtr> cells;
                    CellsGenerator<FixedSimpleWntCellCycleModel, 2> cells_generator;
                    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

                    // Alter cells properties
                    for (unsigned i=0; i<cells.size(); i++)
                    {
                    	// So N(16,1) as in Phil trans paper
                        dynamic_cast<FixedSimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetTransitCellG1Duration(6);
                        dynamic_cast<FixedSimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetWntTransitThreshold(2.0/3.0);
                        dynamic_cast<FixedSimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetWntLabelledThreshold(0.0); // As labelling mutant cells
                    }

                    boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
                    boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

                    // Create cell population
                    PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
                    cell_population.SetNumSweepsPerTimestep(1);

			        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
			        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
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
                    simulator.SetOutputDivisionLocations(true);


                    //Create output directory
                    std::stringstream out;
                    out << index << "/Drag_"<< drag_ratios[drag_index] << "/Height_" << blob_heights[height_index];
                    output_directory = main_directory +  out.str();
                    simulator.SetOutputDirectory(output_directory);

                    simulator.SetDt(0.01);
                    simulator.SetSamplingTimestepMultiple(10);
                    simulator.SetEndTime(time_to_steady_state);
                    simulator.SetOutputCellVelocities(true);

                    // Create cell killer and pass in to simulation
                    MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_length));
                    simulator.AddCellKiller(p_killer);

                    // Create update rules and pass to the simulation
                    MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
                    p_volume_constraint_update_rule->SetMatureCellTargetVolume(25);
                    p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.1); //Default is 0.5
                    simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

                    // DA update rule
                    MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
                    p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.1);
                    p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.2);
                    p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.1);
                    p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.2);
                    p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2);
                    simulator.AddPottsUpdateRule(p_differential_adhesion_update_rule);

                    // Moidifier to track base and top of mutant Patch
                    MAKE_PTR(MutantBaseTrackerModifier<2>, p_base_tracker_modifier);
                    simulator.AddSimulationModifier(p_base_tracker_modifier);

                    // Modifier to track shape of cells
                    MAKE_PTR(CellShapeOutputModifier<2>, p_cell_shape_modifier);
                    simulator.AddSimulationModifier(p_cell_shape_modifier);


                    // Run simulation to steady state
                    simulator.Solve();

                    // Now reset and add mutant cells
                    simulator.SetEndTime(time_to_steady_state + time_after_mutations);

                    // Create some mutant cells
                    c_vector<double, 2> blob_centre;
                    blob_centre[0] = ((double) crypt_width)/ 2.0;
                    blob_centre[1] = blob_heights[height_index];

                    // Iterate over all cells, to define the 'blob'
                    PottsBasedCellPopulation<2>* p_static_cast_cell_population = static_cast<PottsBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation()));

                    for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
                    {
                        double dist_from_blob_centre = 0.0;

                        PottsElement<2>* p_element = p_static_cast_cell_population->GetElementCorrespondingToCell(*cell_iter);

                        for (unsigned index = 0; index<p_element->GetNumNodes(); index++)
                        {
                            dist_from_blob_centre += norm_2(p_element->GetNodeLocation(index)-blob_centre);
                        }
                        dist_from_blob_centre /= p_element->GetNumNodes();

                        if (dist_from_blob_centre < blob_radius)
                        {
                            cell_iter->SetMutationState(p_state);
                            // Also label cells as using Diff adhesion
                            cell_iter->AddCellProperty(p_label);

                        	RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                        	// reset all the cell ages as otherwise the mutant patch will divide straight away as originally old cells
                        	cell_iter->SetBirthTime(time_to_steady_state - 16* p_gen->ranf() );

                        }
                    }

                    // Modify movement of mutant cells with a new update rule
                    MAKE_PTR(MutantCellPottsUpdateRule<2>, p_mutant_cell_update_rule);
                    p_mutant_cell_update_rule->SetMutantCellMovementRatio(drag_ratios[drag_index]);
                    simulator.AddPottsUpdateRule(p_mutant_cell_update_rule);

                    // In order to catch runs with poor motility arameters
                    try
                    {
                    	simulator.Solve();
                    }
                    catch (Exception&)
                    {
                        WARNING("Ignore Run");
                        PRINT_VARIABLE(output_directory);
                    }

                    // Reset singletons as in loop
                    SimulationTime::Destroy();
                    RandomNumberGenerator::Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);
                    WntConcentration<2>::Destroy();
                }
            }
        }
    }
};

#endif /*TESTPOTTSCRYPTMUTANT_HPP_*/