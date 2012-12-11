// -*- C++ -*-
/***************************************************************************
 * blitz/array/iter.h  Basic iterator for arrays.
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
#ifndef BZ_ARRAY_H
 #error <blitz/array/iter.h> must be included via <blitz/array.h>
#endif

#ifndef BZ_ARRAY_ITER_H
#define BZ_ARRAY_ITER_H

#ifdef BZ_HAVE_STL
#include <iterator>
#endif

#if defined(BZ_DEBUG)
#define CheckIteratorValidity(X,Y)                                     \
        BZPRECHECK(data_!=0, X " invalid iterator (empty array)");     \
        BZPRECHECK((data_>=beg_+Y && data_<=end_+Y), ((data_<beg_+Y) ? \
            X " invalid iterator (before beginning of array)" :        \
            X " invalid iterator (past end of array)"));
#else
#define CheckIteratorValidity(X,Y)
#endif

BZ_NAMESPACE(blitz)

template<typename T, int N>
class ConstArrayIterator {
private:
    //  Initialization common to begin,end constructors.
    //
    void Init(const Array<T,N>& array) {
        // Making internal copies of these avoids keeping
        // a pointer to the array and doing indirection.
        lbound_ = array.lbound();
        order_ = array.ordering();

        ubound_(0) = array.ubound(0)+1;
        dataincr_(order_(0)) = array.stride(order_(0));
        for (int i=1,r,s=order_(0);i<N;s=r,++i) {
            r = order_(i);
            ubound_(i) = array.ubound(i)+1;
            dataincr_(r) = array.stride(r)-array.extent(s)*array.stride(s);
        }
#if defined(BZ_DEBUG)
        beg_ = array.data();
        end_ = end_value(array)-1;
        if (beg_>end_)
            std::swap(beg_,end_);
#endif
    }

public:
    ConstArrayIterator() : data_(0) { }

    ConstArrayIterator(const Array<T,N>& array) :
        data_(const_cast<T*>(array.data())) {
        Init(array);
        pos_ = lbound_;
    }

    ConstArrayIterator(const Array<T,N>& array, const int) :
        data_(end_value(array)) {
        Init(array);
        pos_ = array.ubound();
        ++pos_(order_(0));
    }

    bool operator==(const ConstArrayIterator<T,N>& x) const
    { return data_ == x.data_; }

    bool operator!=(const ConstArrayIterator<T,N>& x) const
    { return data_ != x.data_; }

    const T& operator*() const {
        CheckIteratorValidity("Attempted to dereference",0);
        return *data_;
    }

    const T* restrict operator->() const {
        CheckIteratorValidity("Attempted to dereference",0);
        return data_;
    }

    ConstArrayIterator<T,N>& operator++();
    ConstArrayIterator<T,N>& operator--();

    ConstArrayIterator<T,N> operator++(int) {
        ConstArrayIterator<T,N> tmp = *this;
        ++(*this);
        return tmp;
    }

    ConstArrayIterator<T,N> operator--(int) {
        ConstArrayIterator<T,N> tmp = *this;
        --(*this);
        return tmp;
    }

    // get the current position of the Array iterator in index space
    const TinyVector<int,N>& position() const {
        CheckIteratorValidity("Array<T,N>::iterator::position() called on",0);
        return pos_;
    }

private:
    TinyVector<int,N> dataincr_, lbound_, ubound_, order_;

    static T* end_value(const Array<T,N>& array) {
        T* endval = const_cast<T*>(array.data()) +
                                   array.stride(array.ordering(0));
        for (int i=0;i<N;++i)
            endval +=  array.stride(i)*(array.extent(i)-1);
        return endval;
    }

protected:
    TinyVector<int,N> pos_;
    T * restrict data_;
#if defined(BZ_DEBUG)
    const T* restrict beg_;
    const T* restrict end_;
#endif
};


template<typename T, int N>
class ArrayIterator : public ConstArrayIterator<T,N> {
private:
    typedef ConstArrayIterator<T,N> T_base;
    using T_base::data_;

#if defined(BZ_DEBUG)
    using T_base::beg_;
    using T_base::end_;
#endif

public:
    ArrayIterator() { }

    ArrayIterator(Array<T,N>& x) : T_base(x) { }

    ArrayIterator(const Array<T,N>& array, const int): T_base(array,0) { }

    T& operator*() const {
        CheckIteratorValidity("Attempted to dereference",0);
        return *data_;
    }

    T* restrict operator->() const {
        CheckIteratorValidity("Attempted to dereference",0);
        return data_;
    }

    ArrayIterator<T,N>& operator++() {
        T_base::operator++();
        return *this;
    }

    ArrayIterator<T,N> operator++(int) {
        ArrayIterator<T,N> tmp = *this;
        ++(*this);
        return tmp;
    }

    ArrayIterator<T,N>& operator--() {
        T_base::operator--();
        return *this;
    }

    ArrayIterator<T,N> operator--(int) {
        ArrayIterator<T,N> tmp = *this;
        --(*this);
        return tmp;
    }
};

template<typename T, int N>
ConstArrayIterator<T,N>& ConstArrayIterator<T,N>::operator++() {
    CheckIteratorValidity("Attempted to increment",0);

    //   The first loop iteration is peeled as it increases performance.
    //   The same improvement can be obtained by telling the compiler that
    //   the test is likely to be true, but this has too many portability issues
    //   for now.

    // With a compiler peeling loops correctly (or with an effective BZ_LIKELY
    // macro, this could be simply written as:
    //
    // for (int i=0;i<N;++i) {
    //     const int r = order_(i);
    //     data_ += dataincr_[r];
    //     if (BZ_LIKELY(++pos_(r)!=ubound_(r)))
    //         return *this;
    //     pos_(r) = lbound_(r);
    // }


    const int r0 = order_(0);
    data_ += dataincr_[r0];
    if (BZ_LIKELY(++pos_(r0)!=ubound_(r0)))
        return *this;
    pos_(r0) = lbound_(r0);

    for (int i=1;i<N;++i) {
        const int r = order_(i);
        data_ += dataincr_[r];
        if (BZ_LIKELY(++pos_(r)!=ubound_(r)))
            return *this;
        pos_(r) = lbound_(r);
    }

    // At this place the value of data_ should match that of the end iterator.
    // Do the proper correction to achieve that.

    for (int i=1;i<N;++i) {
        const int r = order_(i);
        data_ -= dataincr_[r];
        pos_(r) = ubound_(r)-1;
    }
    pos_(r0) = ubound_(r0);

    return *this;
}

template<typename T, int N>
ConstArrayIterator<T,N>& ConstArrayIterator<T,N>::operator--() {
    CheckIteratorValidity("Attempted to decrement",1);

    //   The first loop iteration is peeled as it increases performance.
    //   The same improvement can be obtained by telling the compiler that
    //   the test is likely to be true, but this has too many portability issues
    //   for now.

    // With a compiler peeling loops correctly (or with an effective BZ_LIKELY
    // macro, this could be simply written as:
    //
    // for (int i=0;i<N;++i) {
    //     const int r = order_(i);
    //     data_ -= dataincr_[r];
    //     if (BZ_LIKELY(pos_(r)--!=lbound_(r)))
    //         return *this;
    //     pos_(r) = ubound_(r)-1;
    // }


    const int r0 = order_(0);
    data_ -= dataincr_[r0];
    if (BZ_LIKELY(pos_(r0)--!=lbound_(r0)))
        return *this;
    pos_(r0) = ubound_(r0)-1;

    for (int i=1;i<N;++i) {
        const int r = order_(i);
        data_ -= dataincr_[r];
        if (BZ_LIKELY(pos_(r)--!=lbound_(r)))
            return *this;
        pos_(r) = ubound_(r)-1;
    }

    // At this place the value of data_ should match that of the end iterator.
    // No correction is needed for operator--

    return *this;
}

BZ_NAMESPACE_END


#ifdef BZ_HAVE_STL
// support for std::iterator_traits
BZ_NAMESPACE(std)

template <typename T, int N>
struct iterator_traits< BZ_BLITZ_SCOPE(ConstArrayIterator)<T,N> > {
    typedef bidirectional_iterator_tag         iterator_category;
    typedef T                                  value_type;
    typedef blitz::diffType                    difference_type;
    typedef const T*                           pointer;
    typedef const T&                           reference;
};

template <typename T, int N>
struct iterator_traits< BZ_BLITZ_SCOPE(ArrayIterator)<T,N> > {
    typedef bidirectional_iterator_tag         iterator_category;
    typedef T                                  value_type;
    typedef blitz::diffType                    difference_type;
    typedef T*                                 pointer;
    typedef T&                                 reference;
};

BZ_NAMESPACE_END

#endif // BZ_HAVE_STL

#endif // BZ_ARRAY_ITER_H

