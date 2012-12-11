// -*- C++ -*-
// $Id$

/*
 * Gamma distribution
 *
 * Source: Ahrens, J.H. and Dieter, U., Generating Gamma variates
 * by a modified rejection technique.  Comm. ACM, 25,1 (Jan. 1982)
 * pp. 47-54.
 *
 * This code has been adapted from RANDLIB.C 1.3, by
 * Barry W. Brown, James Lovato, Kathy Russell, and John Venier.
 * Code was originally by Ahrens and Dieter (see above).
 *
 * Adapter's notes:
 * NEEDS_WORK: more precision for literals.
 * NEEDS_WORK: ideally the normal_ member should be driven from
 * the same IRNG as the Gamma object, in the event that independentState
 * is used.  Not clear how this could be accomplished.
 */

#ifndef BZ_RANDOM_GAMMA
#define BZ_RANDOM_GAMMA

#ifndef BZ_RANDOM_UNIFORM
 #include <random/uniform.h>
#endif

#ifndef BZ_RANDOM_NORMAL
 #include <random/normal.h>
#endif

#ifndef BZ_RANDOM_EXPONENTIAL
 #include <random/exponential.h>
#endif

#ifndef BZ_NUMINQUIRE_H
 #include <blitz/numinquire.h>
#endif

BZ_NAMESPACE(ranlib)

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class Gamma : public UniformOpen<T,IRNG,stateTag>
{
public:
    typedef T T_numtype;

    Gamma()
    {
        setMean(1.0);
    }

  explicit Gamma(unsigned int i) :
    UniformOpen<T,IRNG,stateTag>(i)
  {
    setMean(1.0);
  };

    Gamma(T mean)
    {
        setMean(mean);
    }

  Gamma(T mean, unsigned int i) :
    UniformOpen<T,IRNG,stateTag>(i)
  {
    setMean(mean);
  };

    T random();

    void setMean(T mean)
    {
        BZPRECONDITION(mean >= 1.0);
        a = mean;
    }

protected:
    T ranf()
    {
        return UniformOpen<T,IRNG,stateTag>::random();
    }

    T snorm()
    {
        return normal_.random();
    }

    T sexpo()
    {
        return exponential_.random();
    }

    T fsign(T num, T sign)
    {
        /* Transfers sign of argument sign to argument num */

        if ((sign>0.0L && num<0.0L)||(sign<0.0L && num>0.0L))
            return -num;
        else
            return num;
    }

    NormalUnit<T,IRNG,sharedState> normal_;
    ExponentialUnit<T,IRNG,sharedState> exponential_;

    T a;
};

template<typename T, typename IRNG, typename stateTag>
T Gamma<T,IRNG,stateTag>::random()
{
    /*
     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
     */

static T q1 = 4.166669E-2;
static T q2 = 2.083148E-2;
static T q3 = 8.01191E-3;
static T q4 = 1.44121E-3;
static T q5 = -7.388E-5;
static T q6 = 2.4511E-4;
static T q7 = 2.424E-4;
static T a1 = 0.3333333;
static T a2 = -0.250003;
static T a3 = 0.2000062;
static T a4 = -0.1662921;
static T a5 = 0.1423657;
static T a6 = -0.1367177;
static T a7 = 0.1233795;
static T e1 = 1.0;
static T e2 = 0.4999897;
static T e3 = 0.166829;
static T e4 = 4.07753E-2;
static T e5 = 1.0293E-2;
static T aa = 0.0;
static T aaa = 0.0;
static T sqrt32 = 5.656854249492380195206754896838792314280;

/* JJV added b0 to fix rare and subtle bug */
static T sgamma,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p;

    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = snorm();
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = ranf();
    if(d*u <= t*t*t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463+s+0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-1*s;
S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if(x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if(log(1.0-u) <= q) return sgamma;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = sexpo();
    u = ranf();
    u += (u-1.0);
    t = b+fsign(si*e,u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if(t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
/*
 * JJV modified the code through line 115 to handle large Q case
 */
    if(q < 15.0) goto S95;
/*
 * JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
 * JJV so reformulate test at 110 in terms of one EXP, if not too big
 * JJV 87.49823 is close to the largest real which can be
 * JJV exponentiated (87.49823 = log(1.0E38))
 */
    if((q+e-0.5*t*t) > 87.49823) goto S115;
    if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
    goto S115;
S95:
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
S115:
    x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))

     JJV changed B to B0 (which was added to declarations for this)
     JJV in 120 to END to fix rare and subtle bug.
     JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
     JJV Reasons: the state of AA only serves to tell the A >= 1.0
     JJV case if certain A-dependent constants need to be recalculated.
     JJV The A < 1.0 case (here) no longer changes any of these, and
     JJV the recalculation of B (which used to change with an
     JJV A < 1.0 call) is governed by the state of AAA anyway.
    aa = 0.0;
*/
    b0 = 1.0+0.3678794*a;
S130:
    p = b0*ranf();
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(sexpo() < sgamma) goto S130;
    return sgamma;
S140:
    sgamma = -log((b0-p)/ a);
    if(sexpo() < (1.0-a)*log(sgamma)) goto S130;
    return sgamma;

}

BZ_NAMESPACE_END

#endif // BZ_RANDOM_GAMMA
