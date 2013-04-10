// -*- C++ -*-
// $Id$

/*
 * Generate Beta random deviate
 *
 *   Returns a single random deviate from the beta distribution with
 *   parameters A and B.  The density of the beta is
 *             x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
 *
 * The mean is a/(a+b).
 * The variance is ab/((a+b)^2(a+b+1))
 * The rth moment is (a+r-1)^(r)/(a+b+r-1)^(r)
 *   where a^(r) means a * (a-1) * (a-2) * ... * (a-r+1)
 *
 * Method
 *    R. C. H. Cheng
 *    Generating Beta Variates with Nonintegral Shape Parameters
 *    Communications of the ACM, 21:317-322  (1978)
 *    (Algorithms BB and BC)
 *    http://www.acm.org/pubs/citations/journals/toms/1991-17-1/p98-l_ecuyer/
 *
 * This class has been adapted from RANDLIB.C 1.3, by
 * Barry W. Brown, James Lovato, Kathy Russell, and John Venier.
 *
 * Adapter's note (TV): This code has gone from Pascal to Fortran to C.
 * As a result it is a bit of a mess.  Note also that the algorithms were
 * originally designed for 32-bit float, and so some of the constants
 * below have inadequate precision.  This will not cause problems for
 * casual use, but if you are generating millions of beta variates and
 * rely on some convergence property, you may have want to worry
 * about this.
 *
 * NEEDS_WORK: dig out the original paper and determine these constants
 * to precision adequate for 128-bit float.
 * NEEDS_WORK: turn this into structured code.
 */

#ifndef BZ_RANDOM_BETA
#define BZ_RANDOM_BETA

#ifndef BZ_RANDOM_UNIFORM
 #include <random/uniform.h>
#endif

#ifndef BZ_NUMINQUIRE_H
 #include <blitz/numinquire.h>
#endif

BZ_NAMESPACE(ranlib)

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class Beta : public UniformOpen<T,IRNG,stateTag>
{
public:
    typedef T T_numtype;

    Beta(T a, T b)
    {
      setParameters(a, b);
    }

  Beta(T a, T b, unsigned int i) : UniformOpen<T, IRNG, stateTag>(i)
    {
      setParameters(a, b);
    }

    T random();

    void setParameters(T a, T b)
    {
      aa = a;
      bb = b;
      infnty = 0.3 * blitz::huge(T());
      minlog = 0.085 * blitz::tiny(T());
      expmax = log(infnty);
    }

protected:
    T ranf()
    {
        return UniformOpen<T,IRNG,stateTag>::random();
    }

    T aa, bb;
    T infnty, minlog, expmax;
};

template<typename T, typename IRNG, typename stateTag>
T Beta<T,IRNG,stateTag>::random()
{
/* JJV changed expmax (log(1.0E38)==87.49823), and added minlog */

    // TV: The original code had infnty = 1.0E38, minlog = 1.0E-37.
    // I have changed these to 0.3 * huge and 8.5 * tiny.  For IEEE
    // float this works out to 1.020847E38 and 0.9991702E-37.
    // I changed expmax from 87.49823 to log(infnty), which works out
    // to 87.518866 for float.

    static T olda = minlog;
    static T oldb = minlog;
    static T genbet,a,alpha,b,beta,delta,gamma,k1,k2,r,s,t,u1,u2,v,w,y,z;
    static long qsame;

    // This code determines if the aa and bb parameters are unchanged.
    // If so, some code can be skipped.

    qsame = (olda == aa) && (oldb == bb);

    if (!qsame)
    {
      BZPRECHECK((aa > minlog) && (bb > minlog),
          "Beta RNG: parameters must be >= " << minlog << endl
          << "Parameters are aa=" << aa << " and bb=" << bb << endl);
      olda = aa;
      oldb = bb;
    }

    if (!(min(aa,bb) > 1.0))
      goto S100;
/*
     Alborithm BB
     Initialize
*/
    if(qsame) goto S30;
    a = min(aa,bb);
    b = max(aa,bb);
    alpha = a+b;
    beta = sqrt((alpha-2.0)/(2.0*a*b-alpha));
    gamma = a+1.0/beta;
S30:
S40:
    u1 = ranf();
/*
     Step 1
*/
    u2 = ranf();
    v = beta*log(u1/(1.0-u1));
/* JJV altered this */
    if(v > expmax) goto S55;
/*
 * JJV added checker to see if a*exp(v) will overflow
 * JJV S50 _was_ w = a*exp(v); also note here a > 1.0
 */
    w = exp(v);
    if(w > infnty/a) goto S55;
    w *= a;
    goto S60;
S55:
    w = infnty;
S60:
    z = pow(u1,2.0)*u2;
    r = gamma*v-1.3862944;
    s = a+r-w;
/*
     Step 2
*/
    if(s+2.609438 >= 5.0*z) goto S70;
/*
     Step 3
*/
    t = log(z);
    if(s > t) goto S70;
/*
 *   Step 4
 *
 *    JJV added checker to see if log(alpha/(b+w)) will
 *    JJV overflow.  If so, we count the log as -INF, and
 *    JJV consequently evaluate conditional as true, i.e.
 *    JJV the algorithm rejects the trial and starts over
 *    JJV May not need this here since alpha > 2.0
 */
    if(alpha/(b+w) < minlog) goto S40;
    if(r+alpha*log(alpha/(b+w)) < t) goto S40;
S70:
/*
     Step 5
*/
    if(!(aa == a)) goto S80;
    genbet = w/(b+w);
    goto S90;
S80:
    genbet = b/(b+w);
S90:
    goto S230;
S100:
/*
     Algorithm BC
     Initialize
*/
    if(qsame) goto S110;
    a = max(aa,bb);
    b = min(aa,bb);
    alpha = a+b;
    beta = 1.0/b;
    delta = 1.0+a-b;
    k1 = delta*(1.38889E-2+4.16667E-2*b)/(a*beta-0.777778);
    k2 = 0.25+(0.5+0.25/delta)*b;
S110:
S120:
    u1 = ranf();
/*
     Step 1
*/
    u2 = ranf();
    if(u1 >= 0.5) goto S130;
/*
     Step 2
*/
    y = u1*u2;
    z = u1*y;
    if(0.25*u2+z-y >= k1) goto S120;
    goto S170;
S130:
/*
     Step 3
*/
    z = pow(u1,2.0)*u2;
    if(!(z <= 0.25)) goto S160;
    v = beta*log(u1/(1.0-u1));
/*
 *    JJV instead of checking v > expmax at top, I will check
 *    JJV if a < 1, then check the appropriate values
 */
    if(a > 1.0) goto S135;
/*   JJV a < 1 so it can help out if exp(v) would overflow */
    if(v > expmax) goto S132;
    w = a*exp(v);
    goto S200;
S132:
    w = v + log(a);
    if(w > expmax) goto S140;
    w = exp(w);
    goto S200;
S135:
/*   JJV in this case a > 1 */
    if(v > expmax) goto S140;
    w = exp(v);
    if(w > infnty/a) goto S140;
    w *= a;
    goto S200;
S140:
    w = infnty;
    goto S200;
/*
 * JJV old code
 *    if(!(v > expmax)) goto S140;
 *    w = infnty;
 *    goto S150;
 *S140:
 *    w = a*exp(v);
 *S150:
 *    goto S200;
 */
S160:
    if(z >= k2) goto S120;
S170:
/*
     Step 4
     Step 5
*/
    v = beta*log(u1/(1.0-u1));
/*   JJV same kind of checking as above */
    if(a > 1.0) goto S175;
/* JJV a < 1 so it can help out if exp(v) would overflow */
    if(v > expmax) goto S172;
    w = a*exp(v);
    goto S190;
S172:
    w = v + log(a);
    if(w > expmax) goto S180;
    w = exp(w);
    goto S190;
S175:
/* JJV in this case a > 1.0 */
    if(v > expmax) goto S180;
    w = exp(v);
    if(w > infnty/a) goto S180;
    w *= a;
    goto S190;
S180:
    w = infnty;
/*
 *   JJV old code
 *    if(!(v > expmax)) goto S180;
 *    w = infnty;
 *    goto S190;
 *S180:
 *    w = a*exp(v);
 */
S190:
/*
 * JJV here we also check to see if log overlows; if so, we treat it
 * JJV as -INF, which means condition is true, i.e. restart
 */
    if(alpha/(b+w) < minlog) goto S120;
    if(alpha*(log(alpha/(b+w))+v)-1.3862944 < log(z)) goto S120;
S200:
/*
     Step 6
*/
    if(!(a == aa)) goto S210;
    genbet = w/(b+w);
    goto S220;
S210:
    genbet = b/(b+w);
S230:
S220:
    return genbet;
}

BZ_NAMESPACE_END

#endif // BZ_RANDOM_BETA
