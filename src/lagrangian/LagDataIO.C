//
// LagDataIO.C
//
// Created on 09 Jun 2005
//         by Boyce Griffith (boyce@bigboy.verizon.net).
//
// Last modified: <16.Jun.2005 14:33:42 boyce@mstu1.cims.nyu.edu>
//

#include "LagDataIO.h"

#if 0

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// C++ STDLIB INCLUDES
//
#include <fstream>

// STL INCLUDES
//
#include <set>

// BOOST INCLUDES
//
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/utility/chset.hpp>
using namespace boost::spirit;

// SAMRAI INCLUDES
//
#include "tbox/MPI.h"
#include "tbox/PIO.h"

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifdef DEBUG_NO_INLINE
//#include "LagDataIO.I"
//#endif

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LagDataIO::LagDataIO()
{
    return;
}// LagDataIO

LagDataIO::LagDataIO(
    const LagDataIO& from)
{
    return;
}// LagDataIO

LagDataIO::~LagDataIO()
{
    return;
}// ~LagDataIO

LagDataIO& LagDataIO::operator=(
    const LagDataIO& that)
{
    if (this == &that) return(*this);  // check for self-assignment
    
    return(*this);
}// operator=

#define GETLINE(is,str,line_num)                \
    do                                          \
    {                                           \
        getline(is,str);                        \
        ++line_num;                             \
    }                                           \
    while (0)

// Some parsers for different possible input lines.
namespace
{
    bool parse_uint(
        char const* str,
        int& val)
    {
        // this parser reads in a single integer value, followed by an
        // optional "#"-delimited comment.
        return(parse(str,
                     //  Begin grammar
                     (
                         uint_p[assign_a(val)] >>
                         !comment_p("#")
                         )
                     ,
                     //  End grammar
                     blank_p).full);
    }// parse_uint

    bool parse_uint(
        char const* str,
        int& val0,
        int& val1)
    {
        // this parser reads in a two integer values, followed by an
        // optional "#"-delimited comment.
        return(parse(str,
                     //  Begin grammar
                     (
                         uint_p[assign_a(val0)] >>
                         uint_p[assign_a(val1)] >>
                         !comment_p("#")
                         )
                     ,
                     //  End grammar
                     blank_p).full);
    }// parse_uint

    bool parse_point(
        char const* str,
        double X[NDIM])
    {
        // this parser reads in NDIM real numbers (respectively the
        // coordinates of the marker point) followed by an optional
        // "#"-delimited comment.
        return(parse(str,
                     //  Begin grammar
                     (
#if (NDIM == 2)
                         real_p[assign_a(X[0])] >> real_p[assign_a(X[1])] >>
#endif
#if (NDIM == 3)
                         real_p[assign_a(X[0])] >> real_p[assign_a(X[1])] >> real_p[assign_a(X[2])] >>
#endif
                         !comment_p("#")
                         )
                     ,
                     //  End grammar
                     blank_p).full);
    }// parse_point

    bool parse_point(
        char const* str,
        int& idx,
        double X[NDIM])
    {
        // this parser reads in a single integer value (the Lagrangian
        // index) followed by NDIM real numbers (respectively the
        // coordinates of the marker point), followed by an optional
        // "#"-delimited comment.
        return(parse(str,
                     //  Begin grammar
                     (
                         uint_p[assign_a(idx)] >>
#if (NDIM == 2)
                         real_p[assign_a(X[0])] >> real_p[assign_a(X[1])] >>
#endif
#if (NDIM == 3)
                         real_p[assign_a(X[0])] >> real_p[assign_a(X[1])] >> real_p[assign_a(X[2])] >>
#endif
                         !comment_p("#")
                         )
                     ,
                     //  End grammar
                     blank_p).full);
    }// parse_point
}

void LagDataIO::readFibers(
    const string& file_name,
    int& ngroups,
    blitz::Array<int,1>& nfiber_per_group,
    blitz::Array<int,1>& npoint_per_fiber,
    blitz::Array<blitz::Array<int,2>,1>& point_idxs,
    blitz::Array<blitz::Array<double,3>,1>& point_locs,
    const bool compute_idxs)
{
    static const int mpi_root = 0;
    const int mpi_rank = tbox::MPI::getRank();
    
    // Read in the data file on the root process.
    if (mpi_rank == mpi_root)
    {        
        // Open the input file.
        ifstream is(file_name.c_str());
        string str;
        int line_num = 0;
        
        // Get the number of groups of fibers.
        GETLINE(is,str,line_num);
        if (!parse_uint(str.c_str(), ngroups))
        {
            TBOX_ERROR("HeartModel::readFibers()\n" <<
                       "  error on line " << line_num << " of file " << file_name << "\n" <<
                       "  expected input of the form: n  # number of groups of fibers\n" <<
                       endl);
        }
        
        // Resize the blitz::Array objects to accommodate the number
        // of groups of fibers.
        nfiber_per_group.resize(ngroups);
        npoint_per_fiber.resize(ngroups);
        point_idxs.resize(ngroups);
        point_locs.resize(ngroups);
        
        // For each group, determine the number of fibers in the
        // group, and the number of points per fiber.
        for (int ngr = 0; ngr < ngroups; ++ngr)
        {
            GETLINE(is,str,line_num);
            if (!parse_uint(str.c_str(),
                            nfiber_per_group(ngr),
                            npoint_per_fiber(ngr)))
            {
                TBOX_ERROR("LagDataIO::readFibers()\n" <<
                           "  error on line " << line_num << " of file " << file_name << "\n" <<
                           "  expected input of the form: nfg npf  # num fibers/group, num points/fiber\n" <<
                           endl);
            }
            
            if (nfiber_per_group(ngr) < 1)
            {
                TBOX_ERROR("LagDataIO::readFibers()\n" <<
                           "  error on line " << line_num << " of file " << file_name << "\n" <<
                           "  invalid number of fibers: " << nfiber_per_group(ngr) << "\n" <<
                           "  in group number: " << ngr+1 << "\n" <<
                           endl);
            }

            if (npoint_per_fiber(ngr) < 1)
            {
                TBOX_ERROR("LagDataIO::readFibers()\n" <<
                           "  error on line " << line_num << " of file " << file_name << "\n" <<
                           "  invalid number of points per fiber: " << npoint_per_fiber(ngr) << "\n" <<
                           "  in group number: " << ngr+1 << "\n" <<
                           endl);
            }
        }

        // For each fiber, read in the point locations.
        // 
        // Note: When we read in the index of each fiber point, we
        // keep track of all of the fiber point indices to ensure that
        // each fiber point has a unique Lagrangian index.
        int count = -1;
        set<int> encountered_idxs;
        for (int ngr = 0; ngr < ngroups; ++ngr)
        {
            point_idxs(ngr).resize(nfiber_per_group(ngr),npoint_per_fiber(ngr));
            point_locs(ngr).resize(nfiber_per_group(ngr),npoint_per_fiber(ngr),NDIM);

            for (int nf = 0; nf < nfiber_per_group(ngr); ++nf)
            {
                int idx;
                double X[NDIM];
                for (int np = 0; np < npoint_per_fiber(ngr); ++np)
                {
                    if (compute_idxs)
                    {
                        idx = ++count;
                        if (!parse_point(str.c_str(), X))
                        {
                            TBOX_ERROR("LagDataIO::readFibers()\n" <<
                                       "  error on line " << line_num << " of file " << file_name << "\n" <<
#if (NDIM == 2)
                                       "  expected input of the form: x y\n" <<
#endif
#if (NDIM == 3)
                                       "  expected input of the form: x y z\n" <<
#endif
                                       endl);
                        }
                    }
                    else
                    {
                        if (!parse_point(str.c_str(), idx, X))
                        {
                            TBOX_ERROR("LagDataIO::readFibers()\n" <<
                                       "  error on line " << line_num << " of file " << file_name << "\n" <<
#if (NDIM == 2)
                                       "  expected input of the form: idx x y\n" <<
#endif
#if (NDIM == 3)
                                       "  expected input of the form: idx x y z\n" <<
#endif
                                       endl);
                        }
                        
                        if (encountered_idxs.count(idx) > 0)
                        {                    
                            TBOX_ERROR("LagDataIO::readFibers()\n" <<
                                       "  error on line " << line_num << " of file " << file_name << "\n" <<
                                       "  a fiber point with index " << idx << " has already been encountered." <<
                                       endl);
                        }
                        encountered_idxs.insert(idx);
                    }
                    
                    point_idxs(ngr)(nf,np) = idx;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        point_locs(ngr)(nf,np,d) = X[d];
                    }                    
                }
            }
        }
    }

    // Broadcast the fiber data from the root MPI process to the other
    // MPI processes.
    tbox::MPI::bcast(ngroups, mpi_root);

    nfiber_per_group.resizeAndPreserve(ngroups);
    npoint_per_fiber.resizeAndPreserve(ngroups);
    point_idxs.resizeAndPreserve(ngroups);
    point_locs.resizeAndPreserve(ngroups);
    
    int length;
    for (int ngr = 0; ngr < ngroups; ++ngr)
    {
        tbox::MPI::bcast(nfiber_per_group(ngr), mpi_root);
        tbox::MPI::bcast(npoint_per_fiber(ngr), mpi_root);
        
        point_idxs(ngr).resizeAndPreserve(nfiber_per_group(ngr),npoint_per_fiber(ngr));
        point_locs(ngr).resizeAndPreserve(nfiber_per_group(ngr),npoint_per_fiber(ngr),NDIM);
        
        length = npoint_per_fiber(ngr)*nfiber_per_group(ngr);
        tbox::MPI::bcast(point_idxs(ngr).data(), length, mpi_root);
        
        length = sizeof(double)*NDIM*npoint_per_fiber(ngr)*nfiber_per_group(ngr);
        tbox::MPI::bcast((char*)point_locs(ngr).data(), length, mpi_root);
    }
    
    return;
}// readFibers

void LagDataIO::readMarkers(
    const string& file_name,
    int& nclouds,
    blitz::Array<int,1>& nmarks,
    blitz::Array<blitz::Array<int,1>,1>& marker_idxs,
    blitz::Array<blitz::Array<double,2>,1>& marker_locs,
    const bool compute_idxs)
{
    static const int mpi_root = 0;
    const int mpi_rank = tbox::MPI::getRank();
    
    // Read in the data file on the root process.
    if (mpi_rank == mpi_root)
    {        
        // Open the input file.
        ifstream is(file_name.c_str());
        string str;
        int line_num = 0;
        
        // Get the number of clouds of markers.
        GETLINE(is,str,line_num);
        if (!parse_uint(str.c_str(), nclouds))
        {
            TBOX_ERROR("HeartModel::readMarkers()\n" <<
                       "  error on line " << line_num << " of file " << file_name << "\n" <<
                       "  expected input of the form: n  # number of clouds of markers\n" <<
                       endl);
        }
        
        // Resize the blitz::Array objects to accommodate the number
        // of clouds of markers.
        nmarks.resize(nclouds);
        marker_idxs.resize(nclouds);
        marker_locs.resize(nclouds);
        
        // For each cloud, determine the number of markers per cloud;
        // then read in the marker locations.
        // 
        // Note: When we read in the index of each marker, we keep
        // track of all of the marker indices to ensure that each
        // marker has a unique Lagrangian index.
        int count = -1;
        set<int> encountered_idxs;
        for (int ncl = 0; ncl < nclouds; ++ncl)
        {
            GETLINE(is,str,line_num);

            if (!parse_uint(str.c_str(), nmarks(ncl)))
            {
                TBOX_ERROR("LagDataIO::readMarkers()\n" <<
                           "  error on line " << line_num << " of file " << file_name << "\n" <<
                           "  expected input of the form: n  # number of fluid markers in the cloud\n" <<
                           endl);
            }
            
            if (nmarks(ncl) < 1)
            {
                TBOX_ERROR("LagDataIO::readMarkers()\n" <<
                           "  error on line " << line_num << " of file " << file_name << "\n" <<
                           "  invalid number of points: " << nmarks(ncl) << "\n" <<
                           "  in cloud number: " << ncl+1 << "\n" <<
                           endl);
            }

            marker_idxs(ncl).resize(nmarks(ncl));
            marker_locs(ncl).resize(nmarks(ncl),NDIM);

            int idx;
            double X[NDIM];
            for (int np = 0; np < nmarks(ncl); ++np)
            {
                if (compute_idxs)
                {
                    idx = ++count;
                    if (!parse_point(str.c_str(), X))
                    {
                        TBOX_ERROR("LagDataIO::readMarkers()\n" <<
                                   "  error on line " << line_num << " of file " << file_name << "\n" <<
#if (NDIM == 2)
                                   "  expected input of the form: x y\n" <<
#endif
#if (NDIM == 3)
                                   "  expected input of the form: x y z\n" <<
#endif
                                   endl);
                    }
                }
                else
                {
                    if (!parse_point(str.c_str(), idx, X))
                    {
                        TBOX_ERROR("LagDataIO::readMarkers()\n" <<
                                   "  error on line " << line_num << " of file " << file_name << "\n" <<
#if (NDIM == 2)
                                   "  expected input of the form: idx x y\n" <<
#endif
#if (NDIM == 3)
                                   "  expected input of the form: idx x y z\n" <<
#endif
                                   endl);
                    }
                    
                    if (encountered_idxs.count(idx) > 0)
                    {                    
                        TBOX_ERROR("LagDataIO::readMarkers()\n" <<
                                   "  error on line " << line_num << " of file " << file_name << "\n" <<
                                   "  a marker point with index " << idx << " has already been encountered." <<
                                   endl);
                    }
                    encountered_idxs.insert(idx);
                }
                
                marker_idxs(ncl)(np) = idx;
                for (int d = 0; d < NDIM; ++d)
                {
                    marker_locs(ncl)(np,d) = X[d];
                }
            }
        }
    }

    // Broadcast the marker data from the root MPI process to the
    // other MPI processes.
    tbox::MPI::bcast(nclouds, mpi_root);

    nmarks.resizeAndPreserve(nclouds);
    marker_idxs.resizeAndPreserve(nclouds);
    marker_locs.resizeAndPreserve(nclouds);

    int length;
    for (int ncl = 0; ncl < nclouds; ++ncl)
    {
        tbox::MPI::bcast(nmarks(ncl), mpi_root);

        marker_idxs(ncl).resizeAndPreserve(nmarks(ncl));
        marker_locs(ncl).resizeAndPreserve(nmarks(ncl),NDIM);

        length = nmarks(ncl);
        tbox::MPI::bcast(marker_idxs(ncl).data(), length, mpi_root);

        length = sizeof(double)*NDIM*nmarks(ncl);
        tbox::MPI::bcast((char*)marker_locs(ncl).data(), length, mpi_root);
    }
    
    return;
}// readMarkers

#endif

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
