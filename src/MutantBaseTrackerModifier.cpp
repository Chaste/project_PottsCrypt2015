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

#include "MutantBaseTrackerModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
MutantBaseTrackerModifier<DIM>::MutantBaseTrackerModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
MutantBaseTrackerModifier<DIM>::~MutantBaseTrackerModifier()
{
}

template<unsigned DIM>
void MutantBaseTrackerModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CalculateMutantHeight(rCellPopulation);
}

template<unsigned DIM>
void MutantBaseTrackerModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Create output file
    OutputFileHandler output_file_handler( outputDirectory + "/", false);
    mpMutantBaseResultsFile = output_file_handler.OpenOutputFile("baseofmutations.dat");

    // Write Headers
    *mpMutantBaseResultsFile <<  "time \t mutant region base \t mutant region top \n";

    // Calculate before 1st timestep.
    CalculateMutantHeight(rCellPopulation);
}


template<unsigned DIM>
void MutantBaseTrackerModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Close output file.
    mpMutantBaseResultsFile->close();
}

template<unsigned DIM>
void MutantBaseTrackerModifier<DIM>::CalculateMutantHeight(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double min_mutation_height = DBL_MAX;
    double max_mutation_height = 0.0;

    double height;

    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        if (!cell_iter->GetMutationState()->template IsType<WildTypeCellMutationState>()) // Note the template before the IsType as this is a templated method
        {
            height = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[DIM-1];

            if (height < min_mutation_height)
            {
                min_mutation_height = height;
            }
            if (height > max_mutation_height)
            {
                max_mutation_height = height;
            }
        }
    }
    *mpMutantBaseResultsFile <<  SimulationTime::Instance()->GetTime() << "\t" << min_mutation_height << "\t" << max_mutation_height <<"\n";
}


template<unsigned DIM>
void MutantBaseTrackerModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MutantBaseTrackerModifier<1>;
template class MutantBaseTrackerModifier<2>;
template class MutantBaseTrackerModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MutantBaseTrackerModifier)

