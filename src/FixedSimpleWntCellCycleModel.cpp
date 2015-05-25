/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "FixedSimpleWntCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

FixedSimpleWntCellCycleModel::FixedSimpleWntCellCycleModel()
    : SimpleWntCellCycleModel()
{
}

AbstractCellCycleModel* FixedSimpleWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    FixedSimpleWntCellCycleModel* p_model = new FixedSimpleWntCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetUseCellProliferativeTypeDependentG1Duration(mUseCellProliferativeTypeDependentG1Duration);
    p_model->SetWntStemThreshold(mWntStemThreshold);
    p_model->SetWntTransitThreshold(mWntTransitThreshold);
    p_model->SetWntLabelledThreshold(mWntLabelledThreshold);

    return p_model;
}

void FixedSimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;

    // Set up under what level of Wnt stimulus a cell will divide
    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        wnt_division_threshold = mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
    {
        // should be less than healthy values
        wnt_division_threshold = 0.77*mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
    {
        // less than above value
        wnt_division_threshold = 0.155*mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;

        /*******Begin Change *******/
        // Reset G1Duration in case cell was previously differentiated
        if (mG1Duration == DBL_MAX)
        {
        	SetG1Duration();
        }
        /*******End Change *******/
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = mWntLabelledThreshold;
    }

    double wnt_level = GetWntLevel();
    WntConcentrationType wnt_type = GetWntType();

    // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt_division_threshold
    if (wnt_level >= wnt_division_threshold)
    {
        // For a RADIAL Wnt type, override the cell type to StemCellProliferativeType if the Wnt stimulus exceeds a higher threshold
        if ((wnt_type == RADIAL) && (wnt_level > mWntStemThreshold))
        {
            /*
             * This method is usually called within a CellBasedSimulation, after the CellPopulation
             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
             * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
             * would be creating a new CellPropertyRegistry. In this case the cell proliferative
             * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
             * would be incorrect. We must therefore access the CellProliferativeType via the cell's
             * CellPropertyCollection.
             */
            boost::shared_ptr<AbstractCellProperty> p_stem_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_stem_type);
        }
        else
        {
            boost::shared_ptr<AbstractCellProperty> p_transit_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_transit_type);
        }
    }
    else
    {
        // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void FixedSimpleWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    SimpleWntCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedSimpleWntCellCycleModel)
