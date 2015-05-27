#ifndef TESTPOTTSCRYPTLITERATEPAPER_HPP_
#define TESTPOTTSCRYPTLITERATEPAPER_HPP_

/*
 * = CPM simulation of a healthy crypt =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this test we show how Chaste can be used to simulate a Cellular Potts model of a
 * healthy colon crypt. Full details of the computational model can be found in
 * Osborne (2015) "A Multiscale Model of Colorectal Cancer Using the Cellular Potts Framework".
 *
 * This class was used to produce the data for Figures 2 and 3
 *
 * == Including header files ==
 *
 * EMPTYLINE
 *
 * We begin by including the necessary header files. The first ones are common to all cell_based Chaste simulations
 */

#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"

/* This header includes the Modifier to enable cell shape tracking and can be found in the src folder. */
#include "CellShapeOutputModifier.hpp"

/* The remaining headers are covered in the regular Chaste tutorials */
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
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellIdWriter.hpp"
#include "SloughingCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

/* == Running Simulations ==
 *
 * EMPTYLINE
 *
 * First of all, we define the test class.
 */
class TestPottsCrypt : public AbstractCellBasedTestSuite
{
public:

	/* Although called a test this is the way to run simulations in Chaste */
    void TestSinglePottsCrypt() throw (Exception)
    {
    	/* These variables set up the parameters of the simulation. */
        unsigned crypt_length = 100;
        unsigned crypt_width = 50;
        unsigned element_size = 5;
        double dt = 0.01;// 0.04; 0.02; 0.01, 0.005; 0.0025; 0.00125;
        double end_time = 100; //2200


        /* Create a simple 2D `PottsMesh`, this is the spatial information for the cells. */
        PottsMeshGenerator<2> generator(crypt_width, crypt_width/element_size, element_size, crypt_length +10 , crypt_length/element_size, element_size, 1, 1, 1, true, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();


        /* Create the cells. */
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        CellsGenerator<SimpleWntCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        /* Alter the cells properties as required. */
        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<SimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetTransitCellG1Duration(6);
            dynamic_cast<SimpleWntCellCycleModel*>(cells[i]->GetCellCycleModel())->SetWntTransitThreshold(2.0/3.0);
        }

        /* Create a cell population to associate the cells with the `PottsMesh`. */
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetNumSweepsPerTimestep(1); //Default is 1
        cell_population.SetTemperature(0.1); //Default is 0.1

        /* Specify which statistics to output from the simulation. */
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();


        /* Create an instance of a Wnt concentration, this controls the threshold of cell proliferation. */
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        /* Set up cell-based simulation, and set the relevent parameters. */
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(dt);

        if(dt>1.0)
        {
        	simulator.SetSamplingTimestepMultiple(1);
        }
        else
        {
        	simulator.SetSamplingTimestepMultiple((unsigned)(1.0/dt));
        }
        simulator.SetEndTime(end_time);
        simulator.SetOutputCellVelocities(true);

        simulator.SetOutputDirectory("Potts/HomeostaticCylindricalCrypt");

        /* Create cell killer and pass in to simulation, this specifies how cells die,
         * here they are removed at the top of the crypt. */
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_length));
        simulator.AddCellKiller(p_killer);

        /* Create update rules and pass to the simulation, this determines the terms in the Hamiltonian.
         * First the volume constraint.  */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(25);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.1);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        /* Then the adhesion component, note it is simple to add more components. */
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        /* Add a modifier to track the cell shape. */
        MAKE_PTR(CellShapeOutputModifier<2>, p_cell_shape_modifier);
        simulator.AddSimulationModifier(p_cell_shape_modifier);

        /* Finally run the simulation. */
        simulator.Solve();
    }
    /*
     * With the parameters as above the simulation should take a couple of minutes
	 *
	 * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
	 * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/Potts/HomeostaticCylindricalCrypt/results_from_time_0}}}.
	 * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
	 * java executable. You should also select the axes equal option.
	 *
	 * The data to reproduce Figure 3 can be generated by running this simulation for varying {{{dt}}} and for a longer end time.
	 * The data is in {{{.dat}}} files in the output directory given above.
	 *
   	*/
};

#endif /*TESTPOTTSCRYPTLITERATEPAPER_HPP_*/
