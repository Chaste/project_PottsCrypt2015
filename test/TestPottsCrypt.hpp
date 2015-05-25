#ifndef TESTPOTTSCRYPT_HPP_
#define TESTPOTTSCRYPT_HPP_

#include <cxxtest/TestSuite.h>

#include "../src/CellSurfaceAreasWriter.hpp"
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

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellIdWriter.hpp"

#include "CellShapeOutputModifier.hpp"

#include "SloughingCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
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

    void TestSinglePottsCrypt() throw (Exception)
    {
        unsigned crypt_length = 100; //100
        unsigned crypt_width = 50; //50
        unsigned element_size = 5;


        double dt = 0.01;// 0.04; 0.02; 0.01, 0.005; 0.0025; 0.00125;
        double end_time = 2200; //2200

        // Create a simple 2D PottsMesh
         PottsMeshGenerator<2> generator(crypt_width, crypt_width/element_size, element_size, crypt_length +10 , crypt_length/element_size, element_size, 1, 1, 1, true, true);
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
        cell_population.SetTemperature(0.1); //Default 0.1


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

        std::stringstream out;
		out << "/Dt_" << dt;
		std::string output_directory = "Potts/HomeostaticCylindricalCrypt" +  out.str();
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

        MAKE_PTR(CellShapeOutputModifier<2>, p_cell_shape_modifier);
        simulator.AddSimulationModifier(p_cell_shape_modifier);



        // Run simulation
        simulator.Solve();
    }

};

#endif /*TESTPOTTSCRYPT_HPP_*/
