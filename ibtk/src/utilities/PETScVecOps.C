// Filename: PETScVecOps.C
// Created on 23 Jan 2009 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "PETScVecOps.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/namespaces.h>

// QD INCLUDES
#include <qd/dd_real.h>
#include <qd/qd_real.h>

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
struct greater_abs
    : public std::binary_function<double,double,bool>
{
    inline bool
    operator()(
        const double& x,
        const double& y) const
        {
            return ((x < 0.0 ? -x : x) > (y < 0.0 ? -y : y));
        }
};

struct less_abs
    : public std::binary_function<double,double,bool>
{
    inline bool
    operator()(
        const double& x,
        const double& y) const
        {
            return ((x < 0.0 ? -x : x) < (y < 0.0 ? -y : y));
        }
};

template <class InputIterator,class T>
inline T
compensated_accumulate(
    InputIterator start,
    InputIterator finish,
    T init)
{
    T temp, y, s(init), e(0.0);
    for (InputIterator i = start; i != finish; ++i)
    {
        temp = s;
        y = *i + e;
        s = temp + y;
        e = (temp - s) + y;
    }
    return (s + e);
}// compensated_accumulate
}

PETScVecOps::SortMode      PETScVecOps::s_sort_mode      = NO_SORT;
PETScVecOps::PrecisionMode PETScVecOps::s_precision_mode = DOUBLE;
PETScVecOps::SummationMode PETScVecOps::s_summation_mode = RECURSIVE_SUMMATION;

void
PETScVecOps::setFromDatabase(
    Pointer<Database> db)
{
    if (db.isNull()) return;
    const std::string sort_mode_str = db->getStringWithDefault("sort_mode", "NO_SORT");
    const std::string precision_mode_str = db->getStringWithDefault("precision_mode", "DOUBLE");
    const std::string summation_mode_str = db->getStringWithDefault("summation_mode", "RECURSIVE_SUMMATION");

    if (sort_mode_str == "NO_SORT")
    {
        s_sort_mode = NO_SORT;
    }
    else if (sort_mode_str == "SORT_INCREASING_MAGNITUDE")
    {
        s_sort_mode = SORT_INCREASING_MAGNITUDE;
    }
    else if (sort_mode_str == "SORT_DECREASING_MAGNITUDE")
    {
        s_sort_mode = SORT_DECREASING_MAGNITUDE;
    }
    else
    {
        TBOX_ERROR("PETScVecOps::setFromDatabase():\n"
                   << ":  invalid sort_mode: " << sort_mode_str << ".\n"
                   << "   Choices are: NO_SORT, SORT_INCREASING_MAGNITUDE, SORT_DECREASING_MAGNITUDE.\n");
    }

    if (precision_mode_str == "DOUBLE")
    {
        s_precision_mode = DOUBLE;
    }
    else if (precision_mode_str == "DOUBLE_DOUBLE")
    {
        s_precision_mode = DOUBLE_DOUBLE;
    }
    else if (precision_mode_str == "QUAD_DOUBLE")
    {
        s_precision_mode = QUAD_DOUBLE;
    }
    else
    {
        TBOX_ERROR("PETScVecOps::setFromDatabase():\n"
                   << ":  invalid precision_mode: " << precision_mode_str << ".\n"
                   << "   Choices are: DOUBLE, DOUBLE_DOUBLE, QUAD_DOUBLE.\n");
    }

    if (summation_mode_str == "RECURSIVE_SUMMATION")
    {
        s_summation_mode = RECURSIVE_SUMMATION;
    }
    else if (summation_mode_str == "COMPENSATED_SUMMATION")
    {
        s_summation_mode = COMPENSATED_SUMMATION;
    }
    else
    {
        TBOX_ERROR("PETScVecOps::setFromDatabase():\n"
                   << ":  invalid summation_mode: " << summation_mode_str << ".\n"
                   << "   Choices are: RECURSIVE_SUMMATION, COMPENSATED_SUMMATION.\n");
    }
    return;
}// setFromDatabase

void
PETScVecOps::printClassData(
    std::ostream& os)
{
    os << "PETScVecOps::printClassData():\n";
    os << "  s_sort_mode = ";
    if (s_sort_mode == NO_SORT)
    {
        os << "NO_SORT";
    }
    else if (s_sort_mode == SORT_INCREASING_MAGNITUDE)
    {
        os << "SORT_INCREASING_MAGNITUDE";
    }
    else if (s_sort_mode == SORT_DECREASING_MAGNITUDE)
    {
        os << "SORT_DECREASING_MAGNITUDE";
    }
    else
    {
        os << "UNKNOWN";
    }
    os << "\n";
    os << "  s_precision_mode = ";
    if (s_precision_mode == DOUBLE)
    {
        os << "DOUBLE";
    }
    else if (s_precision_mode == DOUBLE_DOUBLE)
    {
        os << "DOUBLE_DOUBLE";
    }
    else if (s_precision_mode == QUAD_DOUBLE)
    {
        os << "QUAD_DOUBLE";
    }
    else
    {
        os << "UNKNOWN";
    }
    os << "\n";
    os << "  s_summation_mode = ";
    if (s_summation_mode == RECURSIVE_SUMMATION)
    {
        os << "RECURSIVE_SUMMATION";
    }
    else if (s_summation_mode == COMPENSATED_SUMMATION)
    {
        os << "COMPENSATED_SUMMATION";
    }
    else
    {
        os << "UNKNOWN";
    }
    os << "\n";
    return;
}// printClassData

#undef __FUNCT__
#define __FUNCT__ "VecSetValues"
PetscErrorCode
PETScVecOps::VecSetValues(
    Vec xin,
    PetscInt ni,
    const PetscInt ix[],
    const PetscScalar y[],
    InsertMode addv)
{
    PetscFunctionBegin;

    PetscErrorCode ierr;
    PetscMPIInt    rank = xin->stash.rank;
    PetscInt       *owners = xin->map->range,start = owners[rank];
    PetscInt       end = owners[rank+1],i,row;
    PetscScalar    *xx;

#ifdef DEBUG_CHECK_ASSERTIONS
    if (xin->stash.insertmode == INSERT_VALUES && addv == ADD_VALUES)
    {
        TBOX_ERROR("PETScVecOps::VecSetValues():\n"
                   << "  You have already inserted values; you cannot now add.\n");
    }
    else if (xin->stash.insertmode == ADD_VALUES && addv == INSERT_VALUES)
    {
        TBOX_ERROR("PETScVecOps::VecSetValues():\n"
                   << "  You have already added values; you cannot now insert.\n");
    }
#endif
    ierr = VecGetArray(xin,&xx);CHKERRQ(ierr);
    xin->stash.insertmode = addv;

    if (addv == INSERT_VALUES)
    {
        for (i = 0; i < ni; ++i)
        {
            if (xin->stash.ignorenegidx && ix[i] < 0) continue;
#ifdef DEBUG_CHECK_ASSERTIONS
            if (ix[i] < 0)
            {
                TBOX_ERROR("PETScVecOps::VecSetValues():\n"
                           << "  Out of range index value " << ix[i] << " cannot be negative.\n");
            }
#endif
            if ((row = ix[i]) >= start && row < end)
            {
                xx[row-start] = y[i];
            }
            else
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                if (ix[i] >= xin->map->N)
                {
                    TBOX_ERROR("PETScVecOps::VecSetValues():\n"
                               << "  Out of range index value " << ix[i] << " maximum " << xin->map->N << ".\n");
                }
#endif
                ierr = VecStashValue_Private(&xin->stash,row,y[i]);CHKERRQ(ierr);
            }
        }
    }
    else
    {
        for (i = 0; i < ni; ++i)
        {
            if (xin->stash.ignorenegidx && ix[i] < 0) continue;
#ifdef DEBUG_CHECK_ASSERTIONS
            if (ix[i] < 0)
            {
                TBOX_ERROR("PETScVecOps::VecSetValues():\n"
                           << "  Out of range index value " << ix[i] << " cannot be negative.\n");
            }
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
            if (ix[i] > xin->map->N)
            {
                TBOX_ERROR("PETScVecOps::VecSetValues():\n"
                           << "  Out of range index value " << ix[i] << " maximum " << xin->map->N << ".\n");
            }
#endif
            ierr = VecStashValue_Private(&xin->stash,ix[i],y[i]);CHKERRQ(ierr);
        }
    }
    ierr = VecRestoreArray(xin,&xx);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}// VecSetValues

#undef __FUNCT__
#define __FUNCT__ "VecSetValuesBlocked"
PetscErrorCode
PETScVecOps::VecSetValuesBlocked(
    Vec xin,
    PetscInt ni,
    const PetscInt ix[],
    const PetscScalar yin[],
    InsertMode addv)
{
    PetscFunctionBegin;

    PetscMPIInt    rank = xin->stash.rank;
    PetscInt       *owners = xin->map->range,start = owners[rank];
    PetscErrorCode ierr;
    PetscInt       end = owners[rank+1],i,row,bs = xin->map->bs,j;
    PetscScalar    *xx,*y = (PetscScalar*)yin;

    ierr = VecGetArray(xin,&xx);CHKERRQ(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (xin->stash.insertmode == INSERT_VALUES && addv == ADD_VALUES)
    {
        TBOX_ERROR("PETScVecOps::VecSetValuesBlocked():\n"
                   << "  You have already inserted values; you cannot now add.\n");
    }
    else if (xin->stash.insertmode == ADD_VALUES && addv == INSERT_VALUES)
    {
        TBOX_ERROR("PETScVecOps::VecSetValuesBlocked():\n"
                   << "  You have already added values; you cannot now insert.\n");
    }
#endif
    xin->stash.insertmode = addv;

    if (addv == INSERT_VALUES)
    {
        for (i = 0; i < ni; ++i)
        {
            if ((row = bs*ix[i]) >= start && row < end)
            {
                for (j = 0; j < bs; ++j)
                {
                    xx[row-start+j] = y[j];
                }
            }
            else
            {
                if (ix[i] < 0) continue;
#ifdef DEBUG_CHECK_ASSERTIONS
                if (ix[i] >= xin->map->N)
                {
                    TBOX_ERROR("PETScVecOps::VecSetValuesBlocked():\n"
                               << "  Out of range index value " << ix[i] << " max " << xin->map->N << ".\n");
                }
#endif
                ierr = VecStashValuesBlocked_Private(&xin->bstash,ix[i],y);CHKERRQ(ierr);
            }
            y += bs;
        }
    }
    else
    {
        for (i = 0; i < ni; ++i)
        {
            if (ix[i] < 0) continue;
#ifdef DEBUG_CHECK_ASSERTIONS
            if (ix[i] > xin->map->N)
            {
                TBOX_ERROR("PETScVecOps::VecSetValuesBlocked():\n"
                           << "  Out of range index value " << ix[i] << " max " << xin->map->N << ".\n");
            }
#endif
            ierr = VecStashValuesBlocked_Private(&xin->bstash,ix[i],y);CHKERRQ(ierr);
            y += bs;
        }
    }
    ierr = VecRestoreArray(xin,&xx);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}// VecSetValuesBlocked

#undef __FUNCT__
#define __FUNCT__ "VecAssemblyBegin"
PetscErrorCode
PETScVecOps::VecAssemblyBegin(
    Vec xin)
{
    PetscFunctionBegin;

    PetscErrorCode ierr;
    PetscInt       *owners = xin->map->range,*bowners,i,bs,nstash,reallocs;
    PetscMPIInt    size;
    InsertMode     addv;
    MPI_Comm       comm = ((PetscObject)xin)->comm;

    ierr = MPI_Allreduce(&xin->stash.insertmode,&addv,1,MPI_INT,MPI_BOR,comm);CHKERRQ(ierr);
    if (addv == (ADD_VALUES|INSERT_VALUES))
    {
        TBOX_ERROR("PETScVecOps::VecAssemblyBegin():\n"
                   << "  Some processors inserted values while others added.\n");
    }
    xin->stash.insertmode = addv; /* in case this processor had no cache */

    bs = xin->map->bs;
    ierr = MPI_Comm_size(((PetscObject)xin)->comm,&size);CHKERRQ(ierr);
    if (!xin->bstash.bowners && xin->map->bs != -1)
    {
        ierr = PetscMalloc((size+1)*sizeof(PetscInt),&bowners);CHKERRQ(ierr);
        for (i = 0; i < size+1; ++i)
        {
            bowners[i] = owners[i]/bs;
        }
        xin->bstash.bowners = bowners;
    }
    else
    {
        bowners = xin->bstash.bowners;
    }
    ierr = VecStashScatterBegin_Private(&xin->stash,owners);CHKERRQ(ierr);
    ierr = VecStashScatterBegin_Private(&xin->bstash,bowners);CHKERRQ(ierr);
    ierr = VecStashGetInfo_Private(&xin->stash,&nstash,&reallocs);CHKERRQ(ierr);
    ierr = PetscInfo2(xin,"Stash has %D entries, uses %D mallocs.\n",nstash,reallocs);CHKERRQ(ierr);
    ierr = VecStashGetInfo_Private(&xin->bstash,&nstash,&reallocs);CHKERRQ(ierr);
    ierr = PetscInfo2(xin,"Block-Stash has %D entries, uses %D mallocs.\n",nstash,reallocs);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}// VecAssemblyBegin

#undef __FUNCT__
#define __FUNCT__ "VecAssemblyEnd"
PetscErrorCode
PETScVecOps::VecAssemblyEnd(
    Vec xin)
{
    PetscFunctionBegin;

    PetscErrorCode ierr;
    PetscInt       base,i,j,*row,flg,bs;
    PetscMPIInt    rank = xin->stash.rank;
    PetscInt       *owners = xin->map->range,start = owners[rank],end = owners[rank+1];
    PetscMPIInt    n;
    PetscScalar    *val,*vv,*array,*xarray;
    std::vector<std::vector<double> > values(end-start);

    ierr = VecGetArray(xin,&xarray);CHKERRQ(ierr);
    base = xin->map->range[xin->stash.rank];
    bs   = xin->map->bs;

    // Process the stash.
    while (1)
    {
        ierr = VecStashScatterGetMesg_Private(&xin->stash,&n,&row,&val,&flg);CHKERRQ(ierr);
        if (!flg) break;
        if (xin->stash.insertmode == ADD_VALUES)
        {
            for (i = 0; i < n; ++i)
            {
                values[row[i] - base].push_back(val[i]);
            }
        }
        else if (xin->stash.insertmode == INSERT_VALUES)
        {
            for (i = 0; i < n; ++i)
            {
                xarray[row[i] - base] = val[i];
            }
        }
        else
        {
            TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                       << "  Insert mode is not set correctly; corrupted vector.\n");
        }
    }
    ierr = VecStashScatterEnd_Private(&xin->stash);CHKERRQ(ierr);

    // Process the block stash.
    while (1)
    {
        ierr = VecStashScatterGetMesg_Private(&xin->bstash,&n,&row,&val,&flg);CHKERRQ(ierr);
        if (!flg) break;
        for (i = 0; i < n; ++i)
        {
            array = xarray+row[i]*bs-base;
            vv    = val+i*bs;
            if (xin->stash.insertmode == ADD_VALUES)
            {
                for (j = 0; j < bs; ++j)
                {
                    values[row[i]*bs-base+j].push_back(vv[j]);
                }
            }
            else if (xin->stash.insertmode == INSERT_VALUES)
            {
                for (j = 0; j < bs; ++j)
                {
                    array[j] = vv[j];
                }
            }
            else
            {
                TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                           << "  Insert mode is not set correctly; corrupted vector.\n");
            }
        }
    }
    ierr = VecStashScatterEnd_Private(&xin->bstash);CHKERRQ(ierr);

    // Accumulate values.
    if (xin->stash.insertmode == ADD_VALUES)
    {
        for (i = 0; i < end-start; ++i)
        {
            std::vector<double>& i_vals = values[i];
            if (i_vals.empty()) continue;

            // Optionally sort the values.
            if (s_sort_mode == NO_SORT)
            {
                // intentionally blank
            }
            else if (s_sort_mode == SORT_INCREASING_MAGNITUDE)
            {
                std::sort(i_vals.begin(), i_vals.end(), less_abs());
            }
            else if (s_sort_mode == SORT_DECREASING_MAGNITUDE)
            {
                std::sort(i_vals.begin(), i_vals.end(), greater_abs());
            }
            else
            {
                TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                           << "  invalid sort mode; s_sort_mode = " << s_sort_mode << ".\n");
            }

            // Sum the values.
            if (s_precision_mode == DOUBLE)
            {
                if (s_summation_mode == RECURSIVE_SUMMATION)
                {
                    xarray[i] = std::accumulate(i_vals.begin(), i_vals.end(), xarray[i]);
                }
                else if (s_summation_mode == COMPENSATED_SUMMATION)
                {
                    xarray[i] = compensated_accumulate(i_vals.begin(), i_vals.end(), xarray[i]);
                }
                else
                {
                    TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                               << "  invalid summation mode; s_summation_mode = " << s_summation_mode << ".\n");
                }
            }
            else if (s_precision_mode == DOUBLE_DOUBLE)
            {
                unsigned int old_cw;
                fpu_fix_start(&old_cw);
                const std::vector<dd_real> dd_i_vals(i_vals.begin(),i_vals.end());
                if (s_summation_mode == RECURSIVE_SUMMATION)
                {
                    xarray[i] = to_double(std::accumulate(dd_i_vals.begin(), dd_i_vals.end(), dd_real(xarray[i])));
                }
                else if (s_summation_mode == COMPENSATED_SUMMATION)
                {
                    xarray[i] = to_double(compensated_accumulate(dd_i_vals.begin(), dd_i_vals.end(), dd_real(xarray[i])));
                }
                else
                {
                    TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                               << "  invalid summation mode; s_summation_mode = " << s_summation_mode << ".\n");
                }
                fpu_fix_end(&old_cw);
            }
            else if (s_precision_mode == QUAD_DOUBLE)
            {
                unsigned int old_cw;
                fpu_fix_start(&old_cw);
                const std::vector<qd_real> qd_i_vals(i_vals.begin(),i_vals.end());
                if (s_summation_mode == RECURSIVE_SUMMATION)
                {
                    xarray[i] = to_double(std::accumulate(qd_i_vals.begin(), qd_i_vals.end(), qd_real(xarray[i])));
                }
                else if (s_summation_mode == COMPENSATED_SUMMATION)
                {
                    xarray[i] = to_double(compensated_accumulate(qd_i_vals.begin(), qd_i_vals.end(), qd_real(xarray[i])));
                }
                else
                {
                    TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                               << "  invalid summation mode; s_summation_mode = " << s_summation_mode << ".\n");
                }
                fpu_fix_end(&old_cw);
            }
            else
            {
                TBOX_ERROR("PETScVecOps::VecAssemblyEnd():\n"
                           << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
            }
        }
    }

    ierr = VecRestoreArray(xin,&xarray);CHKERRQ(ierr);
    xin->stash.insertmode = NOT_SET_VALUES;

    PetscFunctionReturn(0);
}// VecAssemblyEnd

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
