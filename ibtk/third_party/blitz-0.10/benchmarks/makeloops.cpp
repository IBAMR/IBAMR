

#include <iostream.h>
#include <string.h>
#include <fstream.h>
#include <stdlib.h>
#include <strstream.h>
#include <iostream.h>
#include <stdio.h>
#include <ctype.h>

class loop
{

public:
	loop()
	{
		reset();
	}

	void reset();

	void parseLoop(istream& is);

	int numArrays() const
	{
		return numArrays_;
	}

	char arrayName(int i) const
	{
		return arrays_[i];
	}

	int numScalars() const
	{
		return numScalars_;
	}

	char scalarName(int i) const
	{
		return scalars_[i];
	}

	const char* loopName() const
	{
		return loopName_;
	}

	const char* loopBuffer() const
	{
		return loopBuffer_;
	}

	int flops() const
	{
		return flops_;
	}

	int isArray(char c) const
	{
		for (int i=0; i < numArrays_; ++i)
			if (arrays_[i] == c)
				return 1;
		return 0;
	}

	int isScalar(char c) const
	{
		for (int i=0; i < numScalars_; ++i)
			if (scalars_[i] == c)
				return 1;
		return 0;
	}

private:
	enum { maxArrays = 20, maxScalars = 20, bufLen = 128 };

	char arrays_[maxArrays];
	int numArrays_;
	char scalars_[maxScalars];
	int numScalars_;
	char loopBuffer_[bufLen];
	char loopName_[bufLen];
	int flops_;
};

void loop::reset()
{
	numArrays_ = 0;
	numScalars_ = 0;
	loopBuffer_[0] = '\0';
	loopName_[0] = '\0';
	flops_ = 0;
}

void loop::parseLoop(istream& is)
{
	const int bufLen = 128;
	char buffer[bufLen];

	const char* whitespace = " \t";

	reset();

	while (!is.eof()) {
		is.getline(buffer, bufLen);
		char* token = strtok(buffer, whitespace);
		if (!token)
			continue;

		if (!strcmp(token, "begin")) {
			token = strtok(0, whitespace);
			strcpy(loopName_, token);
			cout << "Creating loop: " << loopName_ << endl;
		} else if (!strcmp(token, "end"))
			return;
		else if (!strcmp(token, "array")) {
			while (token = strtok(0, whitespace)) {
				arrays_[numArrays_++] = token[0];
				cout << "Array: " << token[0] << endl;
			}
		} else if (!strcmp(token, "scalar")) {
			while (token = strtok(0, whitespace)) {
				scalars_[numScalars_++] = token[0];
				cout << "Scalar: " << token[0] << endl;
			}
		} else if (!strcmp(token, "flops")) {
			token = strtok(0, whitespace);
			flops_ = atoi(token);
			cout << "Flops: " << flops_ << endl;
		} else if (!strcmp(token, "loop")) {
			loopBuffer_[0] = '\0';
			while (token = strtok(0, whitespace))
				strcat(loopBuffer_, token);
			cout << "Loop: " << loopBuffer_ << endl;
		}
	}
}

void fortranVersion(loop& lp)
{
	const char* numtype = "REAL*8";

	char filename[128];
	sprintf(filename, "%sf.f", lp.loopName());
	ofstream ofs(filename);

	ofs << "      SUBROUTINE " << lp.loopName() << "_F77(N";
	for (int i=0; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i);
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);
	ofs << ")" << endl
	<< "      INTEGER i, N" << endl
	<< "      " << numtype << " " << lp.arrayName(0) << "(N)";

	for (int i=1; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i) << "(N)";
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);

	ofs << endl << endl
	<< "      DO i=1,N" << endl
	<< "          ";

	const char* loopBuffer = lp.loopBuffer();

	for (int i=0; loopBuffer[i]; ++i) {
		if (loopBuffer[i] == ';')
			ofs << endl << "          ";
		else if (loopBuffer[i] != '$')
			ofs << loopBuffer[i];
		else {
			++i;
			if (lp.isArray(loopBuffer[i]))
				ofs << loopBuffer[i] << "(i)";
			else if (lp.isScalar(loopBuffer[i]))
				ofs << loopBuffer[i];
		}
	}

	ofs << endl
	<< "      END DO" << endl;

	ofs << "      RETURN" << endl
	<< "      END" << endl;

	ofs << endl << endl
	<< "      SUBROUTINE " << lp.loopName() << "_F77Overhead(N";

	for (int i=0; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i);
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);
	ofs << ")" << endl
	<< "      INTEGER i, N" << endl
	<< "      " << numtype << " " << lp.arrayName(0) << "(N)";

	for (int i=1; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i) << "(N)";
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);

	ofs << endl
	<< "      RETURN" << endl
	<< "      END" << endl;
}

void fortran90Version(loop& lp)
{
	const char* numtype = "REAL*8";

	char filename[128];
	sprintf(filename, "%sf90.f90", lp.loopName());
	ofstream ofs(filename);

	ofs << "      SUBROUTINE " << lp.loopName() << "_F90(N";
	for (int i=0; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i);
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);
	ofs << ")" << endl
	<< "      INTEGER i, N" << endl
	<< "      " << numtype << " " << lp.arrayName(0) << "(N)";

	for (int i=1; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i) << "(N)";
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);

	ofs << endl << endl
	<< "      ";

	const char* loopBuffer = lp.loopBuffer();

	for (int i=0; loopBuffer[i]; ++i) {
		if (loopBuffer[i] == ';')
			ofs << endl << "        ";
		else if (loopBuffer[i] != '$')
			ofs << loopBuffer[i];
	}

	ofs << endl
	<< "      RETURN" << endl
	<< "      END" << endl;

	ofs << endl << endl
	<< "      SUBROUTINE " << lp.loopName() << "_F90Overhead(N";
	for (int i=0; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i);
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);
	ofs << ")" << endl
	<< "      INTEGER i, N" << endl
	<< "      " << numtype << " " << lp.arrayName(0) << "(N)";

	for (int i=1; i < lp.numArrays(); ++i)
		ofs << ", " << lp.arrayName(i) << "(N)";
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", " << lp.scalarName(i);

	ofs << endl << endl
	<< "      RETURN" << endl
	<< "      END" << endl;
}

void writeFortranDecl(ofstream& ofs, const char* version, loop& lp,
                      const char* numtype);
void VectorVersion(ostream& os, loop& lp, const char* numtype,
                   const char* scalarArgs, const char* scalarArgs2);
void ArrayVersion(ostream& os, loop& lp, const char* numtype,
                  const char* scalarArgs, const char* scalarArgs2);
void ValarrayVersion(ostream& os, loop& lp, const char* numtype,
                     const char* scalarArgs, const char* scalarArgs2);
void F77Version(ostream& os, loop& lp, const char* numtype,
                const char* scalarArgs, const char* scalarArgs2);
void F90Version(ostream& os, loop& lp, const char* numtype,
                const char* scalarArgs, const char* scalarArgs2);

void cppVersion(loop& lp)
{
	const char* numtype = "double";

	char filename[128];
	sprintf(filename, "%s.cpp", lp.loopName());
	ofstream ofs(filename);

	char capsLoopName[128];
	for (int i=0; i <= strlen(lp.loopName()); ++i)
		capsLoopName[i] = toupper(lp.loopName()[i]);

	ofs << "// Generated code (makeloops.cpp) -- do not edit." << endl << endl
	<< "// In KAI C++ 3.2, restrict causes problems for copy propagation."
	<< endl << "// Temporary fix: disable restrict" << endl << endl
	<< "#define BZ_DISABLE_RESTRICT" << endl << endl
	<<
	"#include <blitz/vector.h>\n"
	"#include <blitz/array.h>\n"
	"#include <blitz/rand-uniform.h>\n"
	"#include <blitz/benchext.h>\n"
	"\n"
	"// Generated: " << __FILE__ << " " << __DATE__ << endl << endl <<
	"#ifdef BZ_HAVE_VALARRAY\n"
	" #define BENCHMARK_VALARRAY\n"
	"#endif\n\n"
	"#ifdef BENCHMARK_VALARRAY\n"
	"#include <valarray>\n"
	"#endif\n"
	"\n"
	"BZ_USING_NAMESPACE(blitz)\n"
	"\n"
	"#if defined(BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES)\n"
	" #define " << lp.loopName() << "_f77 " << lp.loopName() << "_f77_\n"
	" #define " << lp.loopName() << "_f77overhead " << lp.loopName() << "_f77overhead_\n"

	" #define " << lp.loopName() << "_f90 " << lp.loopName() << "_f90_\n"
	" #define " << lp.loopName() << "_f90overhead " << lp.loopName() << "_f90overhead_\n"

	"#elif defined(BZ_FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES)\n"
	" #define " << lp.loopName() << "_f77 " << lp.loopName() << "_f77__\n"
	" #define " << lp.loopName() << "_f77overhead " << lp.loopName() << "_f77overhead__\n"

	" #define " << lp.loopName() << "_f90 " << lp.loopName() << "_f90__\n"
	" #define " << lp.loopName() << "_f90overhead " << lp.loopName() << "_f90overhead__\n"

	"#elif defined(BZ_FORTRAN_SYMBOLS_CAPS)\n"
	" #define " << lp.loopName() << "_f77 " << capsLoopName << "_F77\n"
	" #define " << lp.loopName() << "_f77overhead " << capsLoopName << "_F77OVERHEAD\n"
	" #define " << lp.loopName() << "_f90 " << capsLoopName << "_F90\n"
	" #define " << lp.loopName() << "_f90overhead " << capsLoopName << "_F90OVERHEAD\n"
	"#endif\n"
	"\n"
	"extern \"C\" {" << endl;

	writeFortranDecl(ofs, "_f77", lp, numtype);
	writeFortranDecl(ofs, "_f77overhead", lp, numtype);
	writeFortranDecl(ofs, "_f90", lp, numtype);
	writeFortranDecl(ofs, "_f90overhead", lp, numtype);

	ofs << "}" << endl << endl;

	// Create a string with a list of arguments for the scalars
	ostrstream tmpbuf;
	for (int i=0; i < lp.numScalars(); ++i) {
		tmpbuf << ", " << numtype << " " << lp.scalarName(i);
	}
	tmpbuf << '\0';
	const char* scalarArgs = tmpbuf.str();

	ofs << "void VectorVersion(BenchmarkExt<int>& bench"
	<< scalarArgs << ");" << endl
	<< "void ArrayVersion(BenchmarkExt<int>& bench"
	<< scalarArgs << ");" << endl
	<< "void F77Version(BenchmarkExt<int>& bench"
	<< scalarArgs << ");" << endl
	<< "#ifdef FORTRAN_90" << endl
	<< "void F90Version(BenchmarkExt<int>& bench"
	<< scalarArgs << ");" << endl
	<< "#endif" << endl
	<< "#ifdef BENCHMARK_VALARRAY" << endl
	<< "void ValarrayVersion(BenchmarkExt<int>& bench"
	<< scalarArgs << ");" << endl << "#endif" << endl << endl;

	ofs << "void sink() {}\n\n";

	ofs << "int main()\n"
	"{\n"
	"    int numBenchmarks = 5;\n"
	"#ifndef BENCHMARK_VALARRAY\n"
	"    numBenchmarks--;   // No  valarray\n"
	"#endif\n"
	"#ifndef FORTRAN_90\n"
	"    numBenchmarks--;   // No fortran 90\n"
	"#endif\n"

	"\n"
	"    BenchmarkExt<int> bench(\"" << lp.loopName() << ": "
	<< lp.loopBuffer() << "\", numBenchmarks);\n"
	"\n"
	"    const int numSizes = 23;\n"
	"    bench.setNumParameters(numSizes);\n"
	"    bench.setRateDescription(\"Mflops/s\");\n"
	"\n"
	"    Vector<int> parameters(numSizes);\n"
	"    Vector<long> iters(numSizes);\n"
	"    Vector<double> flops(numSizes);\n"
	"\n"
	"    for (int i=0; i < numSizes; ++i)\n"
	"    {\n"
	"        parameters[i] = (int)pow(10.0, (i+1)/4.0);\n"
	"        iters[i] = 10000000L / parameters[i];\n"
	"        if (iters[i] < 2)\n"
	"            iters[i] = 2;\n"
	"        flops[i] = " << lp.flops() << " * parameters[i];\n"
	"    }\n"
	"\n"
	"    bench.setParameterVector(parameters);\n"
	"    bench.setIterations(iters);\n"
	"    bench.setOpsPerIteration(flops);\n"
	"\n"
	"    bench.beginBenchmarking();" << endl << endl;

	// Create literals
	for (int i=0; i < lp.numScalars(); ++i) {
		ofs << "    " << numtype << " " << lp.scalarName(i)
		<< " = 0.39123982498157938742;" << endl;
	}

	ofs << endl;

	ofs.flush();

	// Create a string with a list of arguments for the scalars
	ostrstream tmpbuf2;
	for (int i=0; i < lp.numScalars(); ++i) {
		tmpbuf2 << ", " << lp.scalarName(i);
	}
	tmpbuf2 << '\0';
	char* scalarArgs2 = tmpbuf2.str();

	ofs << "    VectorVersion(bench" << scalarArgs2 << ");" << endl
	<< "    ArrayVersion(bench" << scalarArgs2 << ");" << endl
	<< "    F77Version(bench" << scalarArgs2 << ");" << endl
	<< "#ifdef FORTRAN_90" << endl
	<< "    F90Version(bench" << scalarArgs2 << ");" << endl
	<< "#endif" << endl
	<< "#ifdef BENCHMARK_VALARRAY" << endl
	<< "    ValarrayVersion(bench" << scalarArgs2 << ");" << endl
	<< "#endif" << endl << endl <<
	"    bench.endBenchmarking();\n"
	"\n"
	"    bench.saveMatlabGraph(\"" << lp.loopName() << ".m\");\n"
	"\n"
	"    return 0;\n"
	"}\n\n"
	"template<class T>\n"
	"void initializeRandomDouble(T data, int numElements, int stride = 1)\n"
	"{\n"
	"    static Random<Uniform> rnd;\n"
	"\n"
	"    for (int i=0; i < numElements; ++i)\n"
	"        data[size_t(i*stride)] = rnd.random();\n"
	"}\n"
	"\n"
	"template<class T>\n"
	"void initializeArray(T& array, int numElements)\n"
	"{\n"
	"    static Random<Uniform> rnd;\n"
	"\n"
	"    for (size_t i=0; i < numElements; ++i)\n"
	"        array[i] = rnd.random();\n"
	"}\n\n";

	VectorVersion(ofs, lp, numtype, scalarArgs, scalarArgs2);
	ArrayVersion(ofs, lp, numtype, scalarArgs, scalarArgs2);
	ValarrayVersion(ofs, lp, numtype, scalarArgs, scalarArgs2);
	F77Version(ofs, lp, numtype, scalarArgs, scalarArgs2);
	F90Version(ofs, lp, numtype, scalarArgs, scalarArgs2);
}

void writeFortranDecl(ofstream& ofs, const char* version, loop& lp,
                      const char* numtype)
{
	ofs << "  void " << lp.loopName() << version
	<< "(const int& N";

	for (int i=0; i < lp.numArrays(); ++i)
		ofs << ", " << numtype << "* " << lp.arrayName(i);
	for (int i=0; i < lp.numScalars(); ++i)
		ofs << ", const " << numtype << "& " << lp.scalarName(i);
	ofs << ");" << endl << endl;
}

void VectorVersion(ostream& os, loop& lp, const char* numtype,
                   const char* scalarArgs, const char* scalarArgs2)
{
	os << "void VectorVersion(BenchmarkExt<int>& bench"
	<< scalarArgs << ")\n"
	<<
	"{\n"
	"    bench.beginImplementation(\"Vector<T>\");\n"
	"\n"
	"    while (!bench.doneImplementationBenchmark())\n"
	"    {\n"
	"        int N = bench.getParameter();\n"
	"        cout << \"Vector<T>: N = \" << N << endl;\n"
	"        cout.flush();\n"
	"\n"
	"        long iters = bench.getIterations();\n"
	"\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        Vector<" << numtype << "> " << lp.arrayName(i)
		<< "(N);" << endl
		<< "        initializeRandomDouble(" << lp.arrayName(i) << ".data(), N);" << endl;
	}

	os << endl <<
	"        bench.start();\n"
	"        for (long i=0; i < iters; ++i)\n"
	"        {\n"
	"            ";

	const char* loopBuffer = lp.loopBuffer();

	for (int i=0; loopBuffer[i]; ++i) {
		if (loopBuffer[i] != '$')
			os << loopBuffer[i];
	}

	os << ";" << endl <<
	"            sink();\n";

	os <<
	"        }\n"
	"        bench.stop();\n\n"
	"        bench.startOverhead();\n"
	"        for (long i=0; i < iters; ++i)\n"
	"            sink();\n"
	"        bench.stopOverhead();\n"

	"    }\n"
	"\n"
	"    bench.endImplementation();\n"
	"}" << endl << endl;
}

void ArrayVersion(ostream& os, loop& lp, const char* numtype,
                  const char* scalarArgs, const char* scalarArgs2)
{
	os << "void ArrayVersion(BenchmarkExt<int>& bench"
	<< scalarArgs << ")\n"
	<<
	"{\n"
	"    bench.beginImplementation(\"Array<T,1>\");\n"
	"\n"
	"    while (!bench.doneImplementationBenchmark())\n"
	"    {\n"
	"        int N = bench.getParameter();\n"
	"        cout << \"Array<T,1>: N = \" << N << endl;\n"
	"        cout.flush();\n"
	"\n"
	"        long iters = bench.getIterations();\n"
	"\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        Array<" << numtype << ", 1> " << lp.arrayName(i)
		<< "(N);" << endl
		<< "        initializeRandomDouble(" << lp.arrayName(i) << ".dataFirst(), N);" << endl;
	}

	os << endl <<
	"        bench.start();\n"
	"        for (long i=0; i < iters; ++i)\n"
	"        {\n"
	"            ";

	const char* loopBuffer = lp.loopBuffer();

	for (int i=0; loopBuffer[i]; ++i) {
		if (loopBuffer[i] != '$')
			os << loopBuffer[i];
	}

	os << ";" << endl <<
	"            sink();\n";

	os <<
	"        }\n"
	"        bench.stop();\n\n"
	"        bench.startOverhead();\n"
	"        for (long i=0; i < iters; ++i)\n"
	"            sink();\n"
	"        bench.stopOverhead();\n";
	os <<
	"    }\n"
	"\n"
	"    bench.endImplementation();\n"
	"}" << endl << endl;
}

void F77Version(ostream& os, loop& lp, const char* numtype,
                const char* scalarArgs, const char* scalarArgs2)
{
	os << "void F77Version(BenchmarkExt<int>& bench"
	<< scalarArgs << ")\n"
	"{\n"
	"    bench.beginImplementation(\"Fortran 77\");\n"
	"\n"
	"    while (!bench.doneImplementationBenchmark())\n"
	"    {\n"
	"        int N = bench.getParameter();\n\n"
	"        cout << \"Fortran 77: N = \" << N << endl;\n"
	"        cout.flush();\n\n"
	"        int iters = bench.getIterations();\n"
	"\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        " << numtype << "* " << lp.arrayName(i)
		<< " = new " << numtype << "[N];" << endl
		<< "        initializeRandomDouble(" << lp.arrayName(i)
		<< ", N);" << endl;
	}

	os << endl <<
	"        bench.start();\n"
	"        for (int iter=0; iter < iters; ++iter)\n"
	"            " << lp.loopName() << "_f77(N";

	for (int i=0; i < lp.numArrays(); ++i)
		os << ", " << lp.arrayName(i);
	os << scalarArgs2 << ");\n"
	"        bench.stop();\n\n"
	"        bench.startOverhead();\n"
	"        for (int iter=0; iter < iters; ++iter)\n"
	"            " << lp.loopName() << "_f77overhead(N";
	for (int i=0; i < lp.numArrays(); ++i)
		os << ", " << lp.arrayName(i);
	os << scalarArgs2 << ");\n";

	os << endl <<
	"        bench.stopOverhead();\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        delete [] " << lp.arrayName(i) << ";" << endl;
	}

	os << "    }\n"
	"\n"
	"    bench.endImplementation();\n"
	"}\n" << endl;

}

void F90Version(ostream& os, loop& lp, const char* numtype,
                const char* scalarArgs, const char* scalarArgs2)
{
	os << "#ifdef FORTRAN_90" << endl
	<< "void F90Version(BenchmarkExt<int>& bench"
	<< scalarArgs << ")\n"
	"{\n"
	"    bench.beginImplementation(\"Fortran 90\");\n"
	"\n"
	"    while (!bench.doneImplementationBenchmark())\n"
	"    {\n"
	"        int N = bench.getParameter();\n\n"
	"        cout << \"Fortran 90: N = \" << N << endl;\n"
	"        cout.flush();\n\n"
	"        int iters = bench.getIterations();\n"
	"\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        " << numtype << "* " << lp.arrayName(i)
		<< " = new " << numtype << "[N];" << endl
		<< "        initializeRandomDouble(" << lp.arrayName(i)
		<< ", N);" << endl;
	}

	os << endl <<
	"        bench.start();\n"
	"        for (int iter=0; iter < iters; ++iter)\n"
	"            " << lp.loopName() << "_f90(N";

	for (int i=0; i < lp.numArrays(); ++i)
		os << ", " << lp.arrayName(i);
	os << scalarArgs2 << ");\n"
	"        bench.stop();\n\n"
	"        bench.startOverhead();\n"
	"        for (int iter=0; iter < iters; ++iter)\n"
	"            " << lp.loopName() << "_f90overhead(N";
	for (int i=0; i < lp.numArrays(); ++i)
		os << ", " << lp.arrayName(i);
	os << scalarArgs2 << ");\n";

	os << endl <<
	"        bench.stopOverhead();\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        delete [] " << lp.arrayName(i) << ";" << endl;
	}

	os << "    }\n"
	"\n"
	"    bench.endImplementation();\n"
	"}\n"
	<< "#endif\n" << endl;

}

void ValarrayVersion(ostream& os, loop& lp, const char* numtype,
                     const char* scalarArgs, const char* scalarArgs2)
{
	os << "#ifdef BENCHMARK_VALARRAY" << endl;
	os << "void ValarrayVersion(BenchmarkExt<int>& bench"
	<< scalarArgs << ")\n"
	<<
	"{\n"
	"    bench.beginImplementation(\"valarray<T>\");\n"
	"\n"
	"    while (!bench.doneImplementationBenchmark())\n"
	"    {\n"
	"        int N = bench.getParameter();\n"
	"        cout << \"valarray<T>: N = \" << N << endl;\n"
	"        cout.flush();\n"
	"\n"
	"        long iters = bench.getIterations();\n"
	"\n";

	for (int i=0; i < lp.numArrays(); ++i) {
		os << "        valarray<" << numtype << "> " << lp.arrayName(i)
		<< "(N);" << endl
		<< "        initializeArray(" << lp.arrayName(i) << ", N);" << endl;
	}

	os << endl <<
	"        bench.start();\n"
	"        for (long i=0; i < iters; ++i)\n"
	"        {\n"
	"            ";

	const char* loopBuffer = lp.loopBuffer();

	for (int i=0; loopBuffer[i]; ++i) {
		if (loopBuffer[i] != '$')
			os << loopBuffer[i];
	}

	os << ";" << endl <<
	"            sink();\n";

	os <<
	"        }\n"
	"        bench.stop();\n\n"
	"        bench.startOverhead();\n"
	"        for (long i=0; i < iters; ++i)\n"
	"            sink();\n"
	"        bench.stopOverhead();\n"

	"    }\n"
	"\n"
	"    bench.endImplementation();\n"
	"}" << endl << endl << "#endif" << endl;
}

int main()
{
	ifstream ifs("loops.data");

	//ofstream ofs("makefile.inc");

	loop lp;

	while (!ifs.eof()) {
		lp.parseLoop(ifs);

		if (ifs.eof())
			break;

		/*
		ofs
#ifdef FORTRAN_90
		<< lp.loopName() << "f90.o:\t" << lp.loopName() << "f90.f"
		<< endl << "\t$(F90) $(F90FLAGS) -c " << lp.loopName() << "f90.f"
		<< endl << endl
#endif
		<< lp.loopName() << ":\t" << lp.loopName() << ".o "
		<< lp.loopName() << "f.o "
#ifdef FORTRAN_90
		<< lp.loopName() << "f90.o"
#endif
		<< endl
		<< "\t$(CXX) $(CXXFLAGS) $(LDFLAGS) -o " << lp.loopName() << " "
		<< lp.loopName() << ".o " << lp.loopName() << "f.o "
#ifdef FORTRAN_90
		<< lp.loopName() << "f90.o "
#endif
		<< "$(LIBS)" << endl << endl;
		*/

		fortranVersion(lp);
//#ifdef FORTRAN_90
		fortran90Version(lp);
//#endif
		cppVersion(lp);
	}
	return 0;
}



