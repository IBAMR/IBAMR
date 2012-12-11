#include "testsuite.h"
#include "blitz/array.h"
using namespace blitz;
using namespace std;

// Test that the reference count works correctly when many threads
// clobber it.

#if !defined(BZ_THREADSAFE) || !defined(_REENTRANT)

int main(int, char**)
{
#ifndef BZ_THREADSAFE
  std::cerr << "must use --enable-threadsafe to test thread safety\n";
#endif
#ifndef _REENTRANT
  std::cerr << "pthreads not enabled, can't create threads (you must use\n--disable-cxx-flags-preset and add -pthread to CXXFLAGS if you are using TBB)\n";
#endif
  return 0;
}
#else

#include <pthread.h>
#include <iostream>
const int N=60000;
const int n_threads=50;

void *run_thread(void *arg)
{
  Array<double,1>* a = 
    reinterpret_cast< Array<double,1>* >(arg);
  // loop and make repeated copies of this array
  for(int i=0; i<N; ++i) {
    Array<double,1> b(*a);
    b+=1;
  }
  pthread_exit(NULL);
}

// Turn off locking by passing an argument (any argument) and watch it
// crash (unless you're using TBB)
int main (int argc, char *argv[])
{
  Array<double,1> a(5);
  a.threadLocal(argc>1);
  a=0;


  pthread_t threads[n_threads];

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  int rc;
  std::cout << "Starting threads\n";
  for(int i=0; i<n_threads; ++i) {
    rc = pthread_create(&threads[i], &attr, run_thread, 
			(void *) &a);
    if(rc) {
      cerr << "Failed at creating thread " << i << endl;
      return 1;
    }
  }

  std::cout << "Joining threads\n";

  // join them
  pthread_attr_destroy(&attr);
  for(int i=0; i<n_threads; i++) {
    rc = pthread_join(threads[i], NULL);
    if (rc) {
      cerr << "Failed at joining thread " << i << endl;
      return 1;
    }
  }
  std::cout << "Done, " << a.numReferences() << " refs left\n";

  BZTEST(a.numReferences()==1);

  // if a wasn't updated, something went wrong.
  BZTEST(all(a>0));

  return 0;
}

#endif
