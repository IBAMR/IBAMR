/***************************************************************************
 * blitz/benchext.cc  Methods for Benchmarking class with external control.
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
#ifndef BZ_BENCHEXT_CC
#define BZ_BENCHEXT_CC

#ifndef BZ_BENCHEXT_H
 #error <blitz/benchext.cc> must be included via <blitz/benchext.h>
#endif

//#include <blitz/vector-et.h>

#ifdef BZ_HAVE_STD
 #include <fstream>
 #include <string>
#else
 #include <fstream.h>
 #include <string.h>
#endif

BZ_NAMESPACE(blitz)

template<typename P_parameter>
BenchmarkExt<P_parameter>::BenchmarkExt(const char* name,
    int numImplementations)
{
    BZPRECONDITION(numImplementations > 0);

    description_ = name;
    numImplementations_ = numImplementations;

    implementationDescriptions_.resize(numImplementations);
    parameterDescription_ = "Vector length";

    // we don't really know this is the default, do we...
    setDependentVariable("flops");

    /*
    // Set up default parameters and iterations
    setNumParameters(19);

    // NEEDS_WORK: once pow(X,Y) is supported, can just say
    // parameters_ = pow(10.0, Range(1,20)/4.0);

    for (unsigned i=0; i < numParameters_; ++i)
      parameters_(i) = static_cast<P_parameter>(BZ_MATHFN_SCOPE(pow)(10.0, (i+1)/4.0));

    iterations_ = 5.0e+5 / parameters_;
    flopsPerIteration_ = parameters_;
    */
    // Set up initial state
    state_ = initializing;
    implementationNumber_ = 0;
}

template<typename P_parameter>
BenchmarkExt<P_parameter>::~BenchmarkExt()
{
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::setNumParameters(int numParameters)
{
    //BZPRECONDITION(state_ == initializing);

    numParameters_ = numParameters;

    parameters_.resize(numParameters_);
    iterations_.resize(numParameters_);
    flopsPerIteration_.resize(numParameters_);

    // Set up timer and Mflops array
    times_.resize(numImplementations_, numParameters_);
    instr_.resize(numImplementations_, numParameters_);
    flops_.resize(numImplementations_, numParameters_);
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::setParameterVector(Array<P_parameter,1> parms)
{
    BZPRECONDITION(state_ == initializing);
    BZPRECONDITION(parms.size() == parameters_.size());

    parameters_ = parms;
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::setParameterDescription(const char* string)
{
    parameterDescription_ = string;
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::setIterations(Array<long,1> iters)
{
    BZPRECONDITION(state_ == initializing);

    iterations_ = iters;
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::setOpsPerIteration(Array<double,1>
    flopsPerIteration)
{
    BZPRECONDITION(flopsPerIteration_.size() == flopsPerIteration.size());

    flopsPerIteration_ = flopsPerIteration;
}

/** Set the dependent variable of the measurements. If the independent
    variable is seconds, we output "G<dvar>/s", if it is cycles, we
    output "<dvar/c". */
template<typename P_parameter>
void BenchmarkExt<P_parameter>::setDependentVariable(const char* dvar)
{
  BZPRECONDITION(Timer::indep_var()!="");

  if(timer_.indep_var()=="s") {
    depvar_ = string("G")+dvar+"/s";
    timerconversion_ = 1./1e9;
  }
  else if(timer_.indep_var()=="c") {
    depvar_ = string(dvar)+"/c";
    timerconversion_ = 1.;
  }
  else {
    depvar_ = std::string(dvar)+"/"+timer_.indep_var();
    timerconversion_ = 1;
    }
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::beginBenchmarking()
{
    BZPRECONDITION(state_ == initializing);
    state_ = benchmarking;
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::beginImplementation(const char* description)
{
    // it would really be better if it worked for as many as we give it
    BZPRECONDITION(implementationNumber_ < numImplementations_);
    BZPRECONDITION(state_ == benchmarking);

    implementationDescriptions_[implementationNumber_] = description;

    state_ = benchmarkingImplementation;
    parameterNumber_ = 0;
}

template<typename P_parameter>
bool BenchmarkExt<P_parameter>::doneImplementationBenchmark() const
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    return parameterNumber_ == numParameters_;
}

template<typename P_parameter>
const string& BenchmarkExt<P_parameter>::currentImplementation() const
{
    BZPRECONDITION(implementationNumber_ < numImplementations_);
    return implementationDescriptions_[implementationNumber_];
}


template<typename P_parameter>
P_parameter BenchmarkExt<P_parameter>::getParameter() const
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    BZPRECONDITION(parameterNumber_ < numParameters_);

    return parameters_(parameterNumber_);
}

template<typename P_parameter>
long BenchmarkExt<P_parameter>::getIterations() const
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    BZPRECONDITION(parameterNumber_ < numParameters_);

    return iterations_(parameterNumber_);
}

template<typename P_parameter>
inline void BenchmarkExt<P_parameter>::start()
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    BZPRECONDITION(parameterNumber_ < numParameters_);
    state_ = running;
    timer_.start();
}

template<typename P_parameter>
inline void BenchmarkExt<P_parameter>::stop()
{
    timer_.stop();
    BZPRECONDITION(state_ == running);
    state_ = benchmarkingImplementation;

    times_(int(implementationNumber_), int(parameterNumber_)) = timer_.elapsed();
    instr_(int(implementationNumber_), int(parameterNumber_)) = timer_.instr();
    flops_(int(implementationNumber_), int(parameterNumber_)) = timer_.flops();
    ++parameterNumber_;
}

template<typename P_parameter>
inline void BenchmarkExt<P_parameter>::startOverhead()
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    BZPRECONDITION(parameterNumber_ > 0);
    BZPRECONDITION(parameterNumber_ <= numParameters_);
    state_ = runningOverhead;
    timer_.start();
}

template<typename P_parameter>
inline void BenchmarkExt<P_parameter>::stopOverhead()
{
    BZPRECONDITION(state_ == runningOverhead);
    timer_.stop();

    cout << "\ttimer overhead: " <<
      1.0*timer_.elapsed()/times_(int(implementationNumber_), int(parameterNumber_-1)) << endl;

    times_(int(implementationNumber_), int(parameterNumber_-1)) -=
      timer_.elapsed();
    instr_(int(implementationNumber_), int(parameterNumber_-1)) -=
      timer_.instr();
    flops_(int(implementationNumber_), int(parameterNumber_-1)) -=
      timer_.flops();

    if(times_(int(implementationNumber_), int(parameterNumber_-1))<0) {
      cerr << "\tError: Timer underflow in benchmark " << implementationDescriptions_[implementationNumber_] << " " << parameters_(parameterNumber_-1) << endl;
      times_(int(implementationNumber_), int(parameterNumber_-1)) = blitz::huge(times_(0,0));
    }
    state_ = benchmarkingImplementation;
}

template<typename P_parameter>
inline void BenchmarkExt<P_parameter>::skip()
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    BZPRECONDITION(parameterNumber_ < numParameters_);
    times_(int(implementationNumber_), int(parameterNumber_)) = blitz::quiet_NaN(double());
    ++parameterNumber_;
    state_ = benchmarkingImplementation;
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::endImplementation()
{
    BZPRECONDITION(state_ == benchmarkingImplementation);
    BZPRECONDITION(parameterNumber_ == numParameters_);

    ++implementationNumber_;

    state_ = benchmarking;
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::endBenchmarking()
{
    BZPRECONDITION(state_ == benchmarking);
    BZPRECONDITION(implementationNumber_ == numImplementations_);

    state_ = done;
}

template<typename P_parameter>
double BenchmarkExt<P_parameter>::getMflops(unsigned implementation,
    unsigned parameterNum) const
{
    BZPRECONDITION(state_ == done);
    BZPRECONDITION(implementation < numImplementations_);
    BZPRECONDITION(parameterNum < numParameters_);
    return 1.0*iterations_(parameterNum) * flopsPerIteration_(parameterNum)
      / times_(int(implementation), int(parameterNum)) * timerconversion_;
}

template<typename P_parameter>
double BenchmarkExt<P_parameter>::getinstrperc(int implementation,
					       int parameterNum) const
{
    BZPRECONDITION(state_ == done);
    BZPRECONDITION(implementation < numImplementations_);
    BZPRECONDITION(parameterNum < numParameters_);
    return 1.0*instr_(implementation,parameterNum)/
      times_(int(implementation), int(parameterNum));
}
template<typename P_parameter>
double BenchmarkExt<P_parameter>::getflopsperc(int implementation,
					       int parameterNum) const
{
    BZPRECONDITION(state_ == done);
    BZPRECONDITION(implementation < numImplementations_);
    BZPRECONDITION(parameterNum < numParameters_);
    return 1.0*flops_(implementation,parameterNum)/
      times_(int(implementation), int(parameterNum));
}

template<typename P_parameter>
void BenchmarkExt<P_parameter>::saveMatlabGraph(const char* filename, const char* graphType) const
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

    // This will be a lot simpler once Matlab-style output formatting
    // of vectors & matrices is finished.

    // ofs << "parm = " << parameters_ << ";" << endl << endl;

    ofs << "parm = [ ";
    unsigned i;
    for (i=0; i < numParameters_; ++i)
      ofs << setprecision(12) << double(parameters_(i)) << " ";
    ofs << "]; " << endl << endl;

    ofs << "Mf = [ ";
    for (i=0; i < numParameters_; ++i)
    {
        for (unsigned j=0; j < numImplementations_; ++j)
        {
            ofs << setprecision(12) << getMflops(j,i) << " ";
        }
        if (i != numParameters_ - 1)
            ofs << ";" << endl;
    }
    ofs << "] ;" << endl << endl;

    ofs << graphType << "(parm,Mf), title('" << description_ << "'), " << endl
        << "    xlabel('" << parameterDescription_ << "'), "
        << "ylabel('" << depvar_ << "')\n"
        << "legend(";

    for (unsigned j=0; j < numImplementations_; ++j)
    {
        ofs << "'" << implementationDescriptions_[j] << "'";
        if (j != numImplementations_ - 1)
            ofs << ", ";
    }

    ofs << ")" << endl;
}


template<typename P_parameter>
void BenchmarkExt<P_parameter>::savePylabGraph(const char* filename, const char* graphType) const
{
    BZPRECONDITION(state_ == done);

    ofstream ofs(filename);

    assert(ofs.good());

    ofs << "# This python file generated automatically by class Benchmark\n"
        << "# of the Blitz++ class library.\n"
	<< "from pylab import *\nfrom numpy import *\n"
	<< "clf()\n";

    ofs.setf(ios::scientific);

    ofs << "legnames=[";
    for (unsigned j=0; j < numImplementations_; ++j)
    {
        ofs << "'" << implementationDescriptions_[j] << "'";
        if (j != numImplementations_ - 1)
            ofs << ", ";
    }

    ofs << "]\n\nparm = array([ ";
    unsigned i;
    for (i=0; i < numParameters_; ++i)
      ofs << setprecision(12) << double(parameters_(i)) << ", ";
    ofs << "])\n\n";

    ofs << "Mf = array([[ ";
    for (i=0; i < numParameters_; ++i)
    {
        if(i>0) ofs << ", [ ";
        for (unsigned j=0; j < numImplementations_; ++j)
        {
	  ofs << setprecision(12) << getMflops(j,i);
	  if(j<numImplementations_-1) ofs << ", ";
        }
	ofs << "]";
    }
    ofs << "])" << endl << endl;
#ifdef BZ_HAVE_LIBPAPI
    // add i/c and flops counters
    ofs << "ic = array([[ ";
    for (i=0; i < numParameters_; ++i)
    {
        if(i>0) ofs << ", [ ";
        for (unsigned j=0; j < numImplementations_; ++j)
        {
	  ofs << setprecision(12) << getinstrperc(j,i);
	  if(j<numImplementations_-1) ofs << ", ";
        }
	ofs << "]";
    }
    ofs << "])" << endl << endl;
    ofs << "fc = array([[ ";
    for (i=0; i < numParameters_; ++i)
    {
        if(i>0) ofs << ", [ ";
        for (unsigned j=0; j < numImplementations_; ++j)
        {
	  ofs << setprecision(12) << getflopsperc(j,i);
	  if(j<numImplementations_-1) ofs << ", ";
        }
	ofs << "]";
    }
    ofs << "])" << endl << endl;
#endif

    ofs << graphType << "(parm,Mf)\ntitle('" << description_ << "')\n"
        << "xlabel('" << parameterDescription_ << "')\n"
        << "ylabel('" << depvar_ << "')\n";

    ofs << "legend(legnames)\n";
}

BZ_NAMESPACE_END

#endif // BZ_BENCHEXT_CC
