#ifndef BZ_POINTERSTACK_H
#define BZ_POINTERSTACK_H

#include <blitz/blitz.h>

BZ_NAMESPACE(blitz)

// helper class ConstPointerStack
template<typename P_numtype, int N_rank>
class ConstPointerStack {
public:
    typedef P_numtype                T_numtype;

    void operator=(const ConstPointerStack<P_numtype,N_rank>& rhs)
    {
        for (int i=0; i<N_rank; ++i)
            stack_[i] = rhs.stack_[i];
    }

    const T_numtype*& operator[](int position)
    {
        return stack_[position];
    }

private:
    const T_numtype *                stack_[N_rank];
};

BZ_NAMESPACE_END

#endif
