/***************************************************************************
 * blitz/bench.cc  Benchmarking class methods.
 *
 * $Id$
 *
 * Copyright (C) 1997-2011 Todd Veldhuizen <tveldhui@acm.org>
 *
 * This file is a part of Blitz.
 *
 * Blitz is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * Blitz is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Blitz.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions:          blitz-devel@lists.sourceforge.net
 * Bugs:                 blitz-support@lists.sourceforge.net
 *
 * For more information, please see the Blitz++ Home Page:
 *    https://sourceforge.net/projects/blitz/
 *
 ****************************************************************************/
#ifndef BZ_BENCH_CC
#define BZ_BENCH_CC

#ifndef BZ_BENCH_H
 #error <blitz/bench.cc> must be included via <blitz/bench.h>
#endif

#ifdef BZ_HAVE_STD
#include <fstream>
#include <string>
#else
 #include <fstream.h>
 #include <string.h>
#endif

BZ_NAMESPACE(blitz)

template<typename P_parameter>
Benchmark<P_parameter>::Benchmark(unsigned numImplementations)
{
    state_ = uninitialized;
    numImplementations_ = numImplementations;
    numStoredImplementations_ = 0;
    implementations_ = new BenchmarkImplementation<P_parameter>* [numImplementations_];
    rates_.resize(numImplementations, numParameterSettings());
    Mflops_.resize(numImplementations, numParameterSettings());
}

template<typename P_parameter>
Benchmark<P_parameter>::~Benchmark()
{
    delete [] implementations_;
}

template<typename P_parameter>
void Benchmark<P_parameter>::addImplementation(
    BenchmarkImplementation<P_parameter> * implementation)
{
    BZPRECONDITION(state_ == uninitialized);
    BZPRECONDITION(numStoredImplementations_ < numImplementations_);

    implementations_[numStoredImplementations_++] = implementation;

    if (numStoredImplementations_ == numImplementations_)
        state_ = initialized;
}

template<typename P_parameter>
void Benchmark<P_parameter>::run(ostream& log)
{
    BZPRECONDITION(state_ == initialized);
    state_ = running;

    Timer t;

    for (unsigned j=0; j < numImplementations_; ++j)
    {
        for (unsigned i=0; i < numParameterSettings(); ++i)
        {
            log  << setw(20) << implementations_[j]->implementationName()
                 << " " << setw(8) << getParameterSetting(i) << "  ";
            log.flush();

            implementations_[j]->initialize(getParameterSetting(i));
            implementations_[j]->tickle();

            unsigned long iterations = getIterationSetting(i);

            t.start();
            implementations_[j]->run(iterations);
            t.stop();
            double tm = t.elapsedSeconds();

            t.start();
            implementations_[j]->runOverhead(iterations);
            t.stop();
            double tmOverhead = t.elapsedSeconds();

            rates_(j,i) = iterations / (tm - tmOverhead);
            Mflops_(j,i) = rates_(j,i)
                * implementations_[j]->flopsPerIteration() / 1.0e+6;

            log << setw(10) << (rates_(j,i)/1.0e+6) << " Mops/s ";

            if (implementations_[j]->flopsPerIteration() != 0)
            {
                log << "[" << setw(7) << Mflops_(j,i) << " Mflops]";
            }

            log << endl;
            log.flush();

            implementations_[j]->done();
        }
    }

    state_ = done;
}

template<typename P_parameter>
double Benchmark<P_parameter>::getMflops(unsigned implementation,
    unsigned setting) const
{
    BZPRECONDITION(state_ == done);
    BZPRECONDITION(implementation < numImplementations_);
    BZPRECONDITION(setting < numParameterSettings());

    return Mflops_(implementation, setting);
}

template<typename P_parameter>
double Benchmark<P_parameter>::getRate(unsigned implementation,
    unsigned setting) const
{
    BZPRECONDITION(state_ == done);
    BZPRECONDITION(implementation < numImplementations_);
    BZPRECONDITION(setting < numParameterSettings());

    return rates_(implementation, setting);
}

template<typename P_parameter>
void Benchmark<P_parameter>::saveMatlabGraph(const char* filename) const
{
    BZPRECONDITION(state_ == done);

    {
      //ugly but saveMatlabGraph is coded into all benchmarks
      std::string pyfn(filename);
      pyfn=pyfn.replace(pyfn.find(".m"),2,std::string(".py"),0,3);
      savePylabGraph(pyfn.c_str());
    }

    ofstream ofs(filename);

    assert(ofs.good());

    ofs << "% This matlab file generated automatically by class Benchmark"
        << endl << "% of the Blitz++ class library." << endl << endl;

    ofs.setf(ios::scientific);

    ofs << "parm = [ ";
    int i;
    for (i=0; i < numParameterSettings(); ++i)
        ofs << setprecision(12) << double(getParameterSetting(i)) << " ";
    ofs << "]; " << endl << endl;

    ofs << "Mf = [ ";
    for (i=0; i < numParameterSettings(); ++i)
    {
        for (int j=0; j < numImplementations_; ++j)
        {
            ofs << setprecision(12) << getMflops(j,i) << " ";
        }
        if (i != numParameterSettings()-1)
            ofs << ";" << endl;
    }
    ofs << "] ;" << endl << endl;

    ofs << "semilogx(parm,Mf), title('" << description() << "'), " << endl
        << "    xlabel('" << parameterDescription() << "'), "
        << "ylabel('Mflops')" << endl
        << "legend(";

    for (int j=0; j < numImplementations_; ++j)
    {
        ofs << "'" << implementations_[j]->implementationName()
            << "'";
        if (j != numImplementations_ - 1)
            ofs << ", ";
    }

    ofs << ")" << endl;
}


template<typename P_parameter>
void Benchmark<P_parameter>::savePylabGraph(const char* filename) const
{
    BZPRECONDITION(state_ == done);

    ofstream ofs(filename);

    assert(ofs.good());

    ofs << "# This python file generated automatically by class Benchmark\n"
        << "# of the Blitz++ class library.\n"
	<< "from pylab import *\nfrom numpy import *\n"
	<< "clf()\n";

    ofs.setf(ios::scientific);

    ofs << "parm = array([ ";
    int i;
    for (i=0; i < numParameterSettings(); ++i)
        ofs << setprecision(12) << double(getParameterSetting(i)) << ", ";
    ofs << "])\n\n";

    ofs << "Mf = array([[ ";
    for (i=0; i < numParameterSettings(); ++i)
    {
        if(i>0) ofs << ", [ ";
        for (int j=0; j < numImplementations_; ++j)
        {
            ofs << setprecision(12) << getMflops(j,i) << ", ";
        }
      ofs << "]";
    }
    ofs << "])" << endl << endl;

    ofs << "semilogx(parm,Mf)\ntitle('" << description() << "')\n"
        << "xlabel('" << parameterDescription() << "')\n"
        << "ylabel('Mflops')\n";
        << "legend([";

    for (int j=0; j < numImplementations_; ++j)
    {
        ofs << "'" << implementations_[j]->implementationName()
            << "'";
        if (j != numImplementations_ - 1)
            ofs << ", ";
    }

    ofs << "])\n";
}

BZ_NAMESPACE_END

#endif // BZ_BENCH_CC
