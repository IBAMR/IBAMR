// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_SAMRAIScopedVectorDuplicate_inl_h
#define included_IBTK_SAMRAIScopedVectorDuplicate_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <CellData.h>
#include <CellVariable.h>
#include <SAMRAIVectorReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
template <typename TYPE>
template <typename INPUT_TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::SAMRAIScopedVectorDuplicate(
    const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, INPUT_TYPE> >& vector,
    const std::string& name)
    : SAMRAIScopedVectorDuplicate(checked_dereference(vector), name)
{
}

template <typename TYPE>
template <typename INPUT_TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::SAMRAIScopedVectorDuplicate(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM, INPUT_TYPE>& vector,
    const std::string& name)
{
    if (std::is_same_v<INPUT_TYPE, TYPE>)
    {
        d_vector = vector.cloneVector(name);
    }
    else
    {
        std::string new_name = (name.empty() ? vector.getName() : name);
        d_vector = new SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>(
            new_name, vector.getPatchHierarchy(), vector.getCoarsestLevelNumber(), vector.getFinestLevelNumber());
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        auto context = var_db->getContext("SAMRAIScopedVectorDuplicate");
        std::unordered_set<int> input_cv_data_ids;
        std::unordered_map<int, int> output_cv_data_id_map;
        for (int i = 0; i < vector.getNumberOfComponents(); ++i)
        {
            // Setup the variable.
            auto input_variable = vector.getComponentVariable(i);
            auto input_data_id = vector.getComponentDescriptorIndex(i);
            auto output_variable = cloneVariable<INPUT_TYPE, TYPE>(input_variable, input_data_id);
            auto output_gcw = getGhostCellWidth<INPUT_TYPE>(input_data_id);
            auto output_data_id = var_db->registerVariableAndContext(output_variable, context, output_gcw);

            // Setup the control volume (if any).
            auto output_cv_data_id = IBTK::invalid_index;
            auto input_cv_data_id = vector.getControlVolumeIndex(i);
            if (input_cv_data_id != IBTK::invalid_index)
            {
                if (output_cv_data_id_map.count(input_cv_data_id))
                {
                    output_cv_data_id = output_cv_data_id_map[input_cv_data_id];
                }
                else
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > input_cv_var;
                    var_db->mapIndexToVariable(input_cv_data_id, input_cv_var);
                    auto output_cv_var = cloneVariable<INPUT_TYPE, TYPE>(input_cv_var, input_cv_data_id);
                    auto output_cv_gcw = getGhostCellWidth<INPUT_TYPE>(input_cv_data_id);
                    output_cv_data_id = var_db->registerVariableAndContext(output_cv_var, context, output_cv_gcw);
                    output_cv_data_id_map[input_cv_data_id] = output_cv_data_id;
                    d_managed_cv_data_ids.insert(output_cv_data_id);
                    allocatePatchData(output_cv_data_id);
                    transformPatchData<INPUT_TYPE, TYPE>(input_cv_data_id, output_cv_data_id);
                }
            }
            d_vector->addComponent(output_variable, output_data_id, output_cv_data_id);
        }
    }

    d_vector->allocateVectorData();
    d_vector->setToScalar(TYPE(0.0), /*interior_only*/ false);
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::operator SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>&()
{
    return *d_vector;
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::operator SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >()
{
    return SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >(&*d_vector, false);
}

template <typename TYPE>
std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > >
SAMRAIScopedVectorDuplicate<TYPE>::getComponentVectors() const
{
    // Setup SAMRAIVectorReal objects to correspond to the individual vector
    // components.
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > > comps;
    for (int comp = 0; comp < d_vector->getNumberOfComponents(); ++comp)
    {
        comps.emplace_back(
            new SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>(d_vector->getName() + "_component_" + std::to_string(comp),
                                                           d_vector->getPatchHierarchy(),
                                                           d_vector->getCoarsestLevelNumber(),
                                                           d_vector->getFinestLevelNumber()));
        comps.back()->addComponent(d_vector->getComponentVariable(comp),
                                   d_vector->getComponentDescriptorIndex(comp),
                                   d_vector->getControlVolumeIndex(comp));
    }
    return comps;
}

template <typename TYPE>
template <typename SRC_TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::transformFromVector(const SAMRAIScopedVectorDuplicate<SRC_TYPE>& src_vector)
{
    transformFromVector(src_vector());
    return;
}

template <typename TYPE>
template <typename SRC_TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::transformFromVector(const SAMRAI::solv::SAMRAIVectorReal<NDIM, SRC_TYPE>& src_vector)
{
    for (int i = 0; i < d_vector->getNumberOfComponents(); ++i)
    {
        transformPatchData<SRC_TYPE, TYPE>(src_vector.getComponentDescriptorIndex(i),
                                           d_vector->getComponentDescriptorIndex(i));
    }
    return;
}

template <typename TYPE>
template <typename DST_TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::transformToVector(SAMRAIScopedVectorDuplicate<DST_TYPE>& dst_vector) const
{
    transformToVector(dst_vector());
    return;
}

template <typename TYPE>
template <typename DST_TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::transformToVector(SAMRAI::solv::SAMRAIVectorReal<NDIM, DST_TYPE>& dst_vector) const
{
    for (int i = 0; i < d_vector->getNumberOfComponents(); ++i)
    {
        transformPatchData<TYPE, DST_TYPE>(d_vector->getComponentDescriptorIndex(i),
                                           dst_vector.getComponentDescriptorIndex(i));
    }
    return;
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::~SAMRAIScopedVectorDuplicate()
{
    if (d_vector)
    {
        for (auto cv_data_id : d_managed_cv_data_ids)
        {
            freePatchData(cv_data_id);
        }
        deallocate_vector_data(*d_vector);
        free_vector_components(*d_vector);
    }
}

template <class TYPE>
template <class SRC_TYPE, class DST_TYPE>
SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >
SAMRAIScopedVectorDuplicate<TYPE>::cloneVariable(SAMRAI::hier::Variable<NDIM>* src_variable,
                                                 const int src_data_idx) const
{
    auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    auto patch_descriptor = var_db->getPatchDescriptor();
    std::string src_var_name = patch_descriptor->mapIndexToName(src_data_idx);
    std::string src_var_id_string(SAMRAI::tbox::Utilities::intToString(src_data_idx, 4));
    std::string new_var_name;
    if (src_var_name.find("-mp_clone_of_id=") == std::string::npos)
    {
        new_var_name = d_vector->getName() + "_" + src_var_name + "-mp_clone_of_id=" + src_var_id_string;
    }
    else
    {
        std::string::size_type last_dash = src_var_name.rfind("=");
        new_var_name = d_vector->getName() + "_" + src_var_name.substr(0, last_dash + 1) + src_var_id_string;
    }
    static unsigned int internal_id = 0;
    new_var_name += "-internal_id=" + SAMRAI::tbox::Utilities::intToString(src_data_idx, internal_id++);

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > dst_variable;
    auto src_variable_cc = dynamic_cast<SAMRAI::pdat::CellVariable<NDIM, SRC_TYPE>*>(src_variable);
    if (src_variable_cc)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM, SRC_TYPE> > src_data_cc_factory =
            patch_descriptor->getPatchDataFactory(src_data_idx);
        dst_variable =
            new SAMRAI::pdat::CellVariable<NDIM, DST_TYPE>(new_var_name, src_data_cc_factory->getDefaultDepth());
    }
    else
    {
        TBOX_ERROR("unsupported component data centering\n");
    }
    return dst_variable;
}

template <class TYPE>
template <class DATA_TYPE>
SAMRAI::hier::IntVector<NDIM>
SAMRAIScopedVectorDuplicate<TYPE>::getGhostCellWidth(const int data_id) const
{
    auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    auto patch_descriptor = var_db->getPatchDescriptor();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM, DATA_TYPE> > data_cc_factory =
        patch_descriptor->getPatchDataFactory(data_id);
    if (data_cc_factory)
    {
        return data_cc_factory->getGhostCellWidth();
    }
    else
    {
        TBOX_ERROR("unsupported component data centering\n");
    }
    return SAMRAI::hier::IntVector<NDIM>(0);
}

template <class TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::allocatePatchData(const int data_id) const
{
    for (int ln = d_vector->getCoarsestLevelNumber(); ln <= d_vector->getFinestLevelNumber(); ++ln)
    {
        d_vector->getPatchHierarchy()->getPatchLevel(ln)->allocatePatchData(data_id);
    }
    return;
}

template <class TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::freePatchData(const int data_id) const
{
    for (int ln = d_vector->getCoarsestLevelNumber(); ln <= d_vector->getFinestLevelNumber(); ++ln)
    {
        d_vector->getPatchHierarchy()->getPatchLevel(ln)->deallocatePatchData(data_id);
    }
    return;
}

template <class TYPE>
template <class SRC_TYPE, class DST_TYPE>
void
SAMRAIScopedVectorDuplicate<TYPE>::transformPatchData(const int src_data_id, const int dst_data_id) const
{
    for (int ln = d_vector->getCoarsestLevelNumber(); ln <= d_vector->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
            d_vector->getPatchHierarchy()->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            auto patch = patch_level->getPatch(p());

            auto convert_array_data = [](const SAMRAI::pdat::ArrayData<NDIM, SRC_TYPE>& src_data,
                                         SAMRAI::pdat::ArrayData<NDIM, DST_TYPE>& dst_data)
            {
                auto src_data_sz = src_data.getBox().size() * src_data.getDepth();
                auto dst_data_sz = dst_data.getBox().size() * dst_data.getDepth();
                TBOX_ASSERT(src_data_sz == dst_data_sz);
                std::transform(src_data.getPointer(),
                               src_data.getPointer() + src_data_sz,
                               dst_data.getPointer(),
                               [](const SRC_TYPE& src) { return static_cast<DST_TYPE>(src); });
            };

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, SRC_TYPE> > src_cc_data =
                patch->getPatchData(src_data_id);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, DST_TYPE> > dst_cc_data =
                patch->getPatchData(dst_data_id);
            if (src_cc_data && dst_cc_data)
            {
                convert_array_data(src_cc_data->getArrayData(), dst_cc_data->getArrayData());
            }
            else
            {
                TBOX_ERROR("unsupported data conversion\n");
            }
        }
    }
}

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorDuplicate_inl_h
