#include "testsuite.h"
#include <blitz/bzconfig.h>
#include <blitz/array.h>
#include <random/uniform.h>

BZ_USING_NAMESPACE(blitz)

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>
#endif

int main()
{
#ifdef BZ_HAVE_BOOST_SERIALIZATION

  stringstream s;

  TinyVector<double, 4> v(1.1,2.2,3.3,4.4);
  TinyMatrix<double, 2, 4> m;
  m=1.1,2.2,3.3,4.4;
  Array<double,1> a(5),b(a);

  {  
    a=2.2, -3.3, 4.4, -5.5, 6.6;
    boost::archive::text_oarchive oa(s);

    // write class instance to archive
    oa << v << m;
    oa << a << b;
  }

  cout << s.str() << endl;

  {
    Array<double,1> aa,bb;
    boost::archive::text_iarchive ia(s);
    TinyVector<double,4> vv;
    TinyMatrix<double, 2, 4> mm;
    ia >> vv >> mm;
    ia >> aa >> bb;

    BZTEST(all(v==vv));
    BZTEST(all(m==mm));
    BZTEST(all(a==aa));
    BZTEST(all(b==bb));
    // make sure aa and bb still refer to the same block
    bb(3)=4711;
    BZTEST(aa(3)==bb(3));

   }
  s.clear();

  // test serialization of random states
  {
    ranlib::Uniform<float> r;
    r.seed(42);

    // save state
    boost::archive::text_oarchive oa(s);
    {
      ranlib::MersenneTwister::T_state state(r.getState());
      oa << state;
    }

    const float n1 = r.random();
    cout << n1 << ".... ";
    const float crap = r.random();

    // restore state
    boost::archive::text_iarchive ia(s);
    {
      ranlib::MersenneTwister::T_state state(r.getState());
      ia >> state;
      r.setState(state);
    }

    const float n2 = r.random();
    cout << n2 << endl;

    BZTEST(n1 == n2);
  }

#else
  cout << "No serialization support enabled, can't test\n";
#endif
  
  return 0;
}

