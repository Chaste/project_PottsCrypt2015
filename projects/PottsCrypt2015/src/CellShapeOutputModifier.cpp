/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellShapeOutputModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CellShapeOutputModifier<DIM>::CellShapeOutputModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellShapeOutputModifier<DIM>::~CellShapeOutputModifier()
{
}

template<unsigned DIM>
void CellShapeOutputModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CalculateFractionalLength(rCellPopulation);
}

template<unsigned DIM>
void CellShapeOutputModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Create output file
	OutputFileHandler output_file_handler( outputDirectory + "/", false);
    mpCellShapeResultsFile = output_file_handler.OpenOutputFile("cellshapes.dat");

    // Write Headers
    *mpCellShapeResultsFile <<  "time \t fractional length \t total length \n";

    // Calculate before 1st timestep.
    CalculateFractionalLength(rCellPopulation);
}


template<unsigned DIM>
void CellShapeOutputModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Close output file.
    mpCellShapeResultsFile->close();
}

template<unsigned DIM>
void CellShapeOutputModifier<DIM>::CalculateCellShape(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    *mpCellShapeResultsFile <<  SimulationTime::Instance()->GetTime() << "\t";

    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		PottsMesh<DIM>* p_potts_mesh = (&(dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh()));

		// Loop over cells and find associated elements so in the same order as the cells in output files
		for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
			 cell_iter != pCellPopulation->End();
			 ++cell_iter)
		{

			unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(cell_iter);
			unsigned cell_id = cell_iter->GetCellId();
			c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(cell_iter);
			double surface_area = p_potts_mesh->GetSurfaceAreaOfElement(location_index);
			double volume = p_potts_mesh->GetVolumeOfElement(location_index);

			*mpCellShapeResultsFile << location_index << "\t" << cell_id << "\t";
			for (unsigned i=0; i<SPACE_DIM; i++)
			{
				*mpCellShapeResultsFile << centre_location[i] << "\t";
			}
			*mpCellShapeResultsFile << volume <<surface_area << "\t";
		}
	}


    else
    {
        EXCEPTION("CellShapeOutputModifier only works for Potts simulations at present");
    }

    *mpCellShapeResultsFile << "\n";
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellShapeOutputModifier<1>;
template class CellShapeOutputModifier<2>;
template class CellShapeOutputModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellShapeOutputModifier)

