// -*- C++ -*-
/***************************************************************************
 * blitz/memblock.h      MemoryBlock<T> and MemoryBlockReference<T>
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
 ***************************************************************************/

#ifndef BZ_MEMBLOCK_H
#define BZ_MEMBLOCK_H

#include <blitz/blitz.h>
#include <blitz/simdtypes.h>
#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/collection_size_type.hpp>
#include <boost/serialization/nvp.hpp>
#endif
#include <stddef.h>     // diffType

BZ_NAMESPACE(blitz)

enum preexistingMemoryPolicy {
  duplicateData,
  deleteDataWhenDone,
  neverDeleteData
};

// This function makes it possible for users to make sure threadsafety
// is enabled.
inline bool isThreadsafe() {
#ifdef BZ_THREADSAFE
  return true;
#else
  return false;
#endif
}


// Forward declaration of MemoryBlockReference
template<typename T_type> class MemoryBlockReference;

// Class MemoryBlock provides a reference-counted block of memory.  This block
// may be referred to by multiple vector, matrix and array objects.  The memory
// is automatically deallocated when the last referring object is destructed.
// MemoryBlock may be subclassed to provide special allocators.
template<typename P_type>
class MemoryBlock {

    friend class MemoryBlockReference<P_type>;

public:
    typedef P_type T_type;

protected:
    explicit MemoryBlock(sizeType items)
    {
      // pad the length to vecWidth, if not already done
    const int w = simdTypes<T_type>::vecWidth;
    const int mod = items%w;
    if (mod>0)
      items += simdTypes<T_type>::vecWidth-mod;
    BZASSERT(items%w==0);

        length_ = items;
        allocate(length_);

        BZASSERT(dataBlockAddress_ != 0);

        references_ = 1;
        BZ_MUTEX_INIT(mutex)
    }

  /** Constructor for a preallocated block that should be deleted when
      we are done? */
    MemoryBlock(sizeType length, T_type* data)
    {
      length_ = length;
      data_ = data;
      dataBlockAddress_ = data;
      references_ = 1;
      BZ_MUTEX_INIT(mutex);
      allocatedByUs_ = false;
    }

    virtual ~MemoryBlock()
    {
        if (dataBlockAddress_)
        {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
    cout << "MemoryBlock:     freed " << setw(8) << length_
         << " at " << ((void *)dataBlockAddress_) << endl;
#endif

            deallocate();
        }

        BZ_MUTEX_DESTROY(mutex)
    }

    // set mutex locking policy and return true if successful
    bool doLock(bool lockingPolicy)
    {
#if defined(BZ_THREADSAFE) && !defined(BZ_THREADSAFE_USE_ATOMIC)
        if (mutexLocking_ == lockingPolicy) { // already set
            return true;
        }
        else if (references_ <= 1) { // no multiple references, safe to change
            mutexLocking_ = lockingPolicy;
            return true;
        }
        return false; // unsafe to change
#elif defined(BZ_THREADSAFE_USE_ATOMIC)
	// with locking we consider the request successful
	return true;
#else
	// without threadsafety we return false if user wants to enable locking
	return !lockingPolicy;
#endif
    }

    /* Note that the the MemoryBlock will be created with reference
       count 1, so there is no need to call addReference immediately
       upon creating it. (The creator obviously does have to call
       removeReference, though.) This avoids the initial mutex lock. */
    void          addReference()
    {
#ifdef BZ_DEBUG_LOG_REFERENCES
        BZ_MUTEX_LOCK(mutex)
        const int refcount = ++references_;
        BZ_MUTEX_UNLOCK(mutex)

        cout << "MemoryBlock:    reffed " << setw(8) << length_
	       << " at " << ((void *)dataBlockAddress_) << " (r="
	       << refcount << ")" << endl;
#else
        BZ_MUTEX_LOCK(mutex)
        ++references_;
        BZ_MUTEX_UNLOCK(mutex)
#endif
    }

    T_type* restrict      data()
    {
        return data_;
    }

    const T_type* restrict data()      const
    {
        return data_;
    }

    T_type*&      dataBlockAddress()
    {
        return dataBlockAddress_;
    }

    sizeType        length()    const
    {
        return length_;
    }

    int           removeReference()
    {

        BZ_MUTEX_LOCK(mutex)
        const int refcount = --references_;
        BZ_MUTEX_UNLOCK(mutex)

#ifdef BZ_DEBUG_LOG_REFERENCES
	  cout << "MemoryBlock: dereffed  " << setw(8) << length_
	       << " at " << ((void *)dataBlockAddress_) << " (r=" << refcount
	       << ")" << endl;
#endif
        return refcount;
    }

    int references() const
  {
        BZ_MUTEX_LOCK(mutex)
        const int refcount = references_;
	BZ_MUTEX_UNLOCK(mutex)

        return refcount;
    }

protected:
    inline void allocate(sizeType length);
    void deallocate();

private:   // Disabled member functions
  MemoryBlock(const MemoryBlock<T_type>&) { };

  void operator=(const MemoryBlock<T_type>&) {};

  /** The default constructor is needed for serialization. */
  MemoryBlock() {};
#ifdef BZ_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;


  /** Implementation of boost::serialization. It doesn't make sense to
      save the mutex, the reference count, or the allocation flags,
      since upon restoring, we must necessarily reallocate. */
    template<class T_arch>
    void save(T_arch& ar, const unsigned int version) const {
#if defined(BZ_THREADSAFE) && !defined(BZ_THREADSAFE_USE_ATOMIC)
      ar << mutexLocking_;
#endif
      const boost::serialization::collection_size_type count(length_);
      ar << BOOST_SERIALIZATION_NVP(count);
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
      cout << "MemoryBlock: serializing " << length_ << " data for MemoryBlock at "
	   << ((void*)this) << endl;
#endif
      if(length_>0)
	ar << boost::serialization::make_array(data_, length_);
    };

    template<class T_arch>
    void load(T_arch& ar, const unsigned int version) {
#if defined(BZ_THREADSAFE) && !defined(BZ_THREADSAFE_USE_ATOMIC)
      ar >> mutexLocking_;
#endif
      boost::serialization::collection_size_type count(length_);
      ar >> BOOST_SERIALIZATION_NVP(count);
      length_=count;
      allocate(length_);

#ifdef BZ_DEBUG_LOG_ALLOCATIONS
      cout << "MemoryBlock: unserializing " << length_ << " data for MemoryBlock at "
	   << ((void*)this) << endl;
#endif

      if(length_>0)
	ar >> boost::serialization::make_array(data_, length_);

      // initialize the members that are not restored. Note that we
      // must initialize references_ to 0 here because the
      // unserialization always adds a reference since the creation of
      // the MemoryBlock object happens by the boost::serialization
      // magic.
      references_ = 0;
      BZ_MUTEX_INIT(mutex);
    };

  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif


private:   // Data members
#if defined(BZ_THREADSAFE) && !defined(BZ_THREADSAFE_USE_ATOMIC)
    // with atomic reference counts, there is no locking
    bool    mutexLocking_;
#endif
  /** Keeps track of whether the block was preallocated or not. This
      affects potential alignment so must be taken into account when
      we delete. */
  bool allocatedByUs_;

 union {
    T_type * restrict     data_;
   typename simdTypes<T_type>::vecType * restrict data_tv_;
   char * restrict        data_char_;
  };
 union {
   T_type *              dataBlockAddress_;
   typename simdTypes<T_type>::vecType * restrict dBA_tv_;
   char *                dBA_char_;
 };
  sizeType              length_;

    BZ_REFCOUNT_DECLARE(references_)
    BZ_MUTEX_DECLARE(mutex)
};




template<typename P_type>
class MemoryBlockReference {

public:
    typedef P_type T_type;

protected:
    T_type * restrict data_;

private:
    MemoryBlock<T_type>* block_;

#ifdef BZ_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;

  /** Serialization operator. This is a bit hacky, because we need to
      restore the data pointer as part of the skeleton, not
      content. For this reason, we serialize the offset as a
      collection_size_item even though it's a signed type. This makes
      the serialization code treat it as part of the skeleton and not
      content, which is what we want.
  */
    template<class T_arch>
    void save(T_arch& ar, const unsigned int version) const {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
      cout << "MemoryBlockReference: serializing object at " << ((void*)this) << ", MemoryBlock at "
	     << ((void*)block_) <<endl;
#endif
      ar << block_;
      ptrdiff_t ptroffset=0;
      if(block_)
	ptroffset = data_ - block_->data();
      boost::serialization::collection_size_type
	offset(*reinterpret_cast<size_t*>(&ptroffset));
      ar << BOOST_SERIALIZATION_NVP(offset);
    };

  /** Unserialization operator. See comment for the save method for
      the reinterpret_cast hack.
  */
    template<class T_arch>
    void load(T_arch& ar, const unsigned int version) {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
      cout << "MemoryBlockReference: unserializing at " << ((void*)this) << endl;
#endif
      ar >> block_;
      addReference();
      boost::serialization::collection_size_type offset;
      ar >> BOOST_SERIALIZATION_NVP(offset);
      ptrdiff_t ptroffset = *reinterpret_cast<ptrdiff_t*>(&offset);
      if(block_)
	data_ = block_->data() + ptroffset;
      else
	data_ = 0;
    };

  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

public:

    MemoryBlockReference()
    {
        block_ = 0;
	// no block, so nothing to add reference to
        data_ = 0;
    }

    MemoryBlockReference(MemoryBlockReference<T_type>& ref, sizeType offset=0)
    {
        block_ = ref.block_;
        addReference();
        data_ = ref.data_ + offset;
    }

    MemoryBlockReference(sizeType length, T_type* data,
        preexistingMemoryPolicy deletionPolicy)
    {
        // Create a memory block using already allocated memory.

        // Note: if the deletionPolicy is duplicateData, this must
        // be handled by the leaf class.  In MemoryBlockReference,
        // this is treated as neverDeleteData; the leaf class (e.g. Array)
        // must duplicate the data.

        if ((deletionPolicy == neverDeleteData)
            || (deletionPolicy == duplicateData)) {
	    // in this case, we do not need a MemoryBlock to ref-count the data
            block_ = 0;
        }
        else if (deletionPolicy == deleteDataWhenDone) {
            block_ = new MemoryBlock<T_type>(length, data);

#ifdef BZ_DEBUG_LOG_ALLOCATIONS
            cout << "MemoryBlockReference: created MemoryBlock at "
                 << ((void*)block_) << endl;
#endif
        }
	// creating a MemoryBlock automatically sets it to one
	// reference, so we do no longer need to add a reference in
	// the constructor.

        data_ = data;
    }

    explicit MemoryBlockReference(sizeType items)
    {
        block_ = new MemoryBlock<T_type>(items);
	// creating a MemoryBlock automatically sets it to one
	// reference, so we do no longer need to add a reference in
	// the constructor.
        data_ = block_->data();

#ifdef BZ_DEBUG_LOG_ALLOCATIONS
	cout << "MemoryBlockReference: created MemoryBlock at "
	     << ((void*)block_) << endl;
#endif
    }

   ~MemoryBlockReference()
    {
        blockRemoveReference();
    }

  /** Returns true if the offset from data_ is vector aligned. */
  bool isVectorAligned(size_t offset) const
  { return simdTypes<T_type>::isVectorAligned(data_ + offset); }

  /** Returns the allocated length of the memory block. */
  sizeType blockLength() const { return block_->length(); };

protected:
#ifdef BZ_TESTSUITE
public:
#endif
    int numReferences() const
    {
        if (block_)
            return block_->references();
#ifdef BZ_DEBUG_LOG_REFERENCES
        cout << "Invalid reference count for data at "<< data_ << endl;
#endif
        return -1;
    }

    bool lockReferenceCount(bool lockingPolicy) const
    {
        if (block_)
            return block_->doLock(lockingPolicy);
        // if we have no block, consider request successful
#ifdef BZ_DEBUG_LOG_REFERENCES
        cout << "No reference count locking for data at "<< data_ << endl;
#endif
        return true;
    }

    void changeToNullBlock()
    {
        blockRemoveReference();
        block_ = 0;
	// no block, so nothing to add reference to
        data_ = 0;
    }

    void changeBlock(MemoryBlockReference<T_type>& ref, sizeType offset=0)
    {
        blockRemoveReference();
        block_ = ref.block_;
        addReference();
        data_ = ref.data_ + offset;
    }

    void newBlock(sizeType items)
    {
        blockRemoveReference();
        block_ = new MemoryBlock<T_type>(items);
	// creating a memory block automatically sets it to one reference
        data_ = block_->data();

#ifdef BZ_DEBUG_LOG_ALLOCATIONS
	cout << "MemoryBlockReference: created MemoryBlock at "
	     << ((void*)block_) << endl;
#endif
    }

private:
    void blockRemoveReference()
    {
        const int refcount = removeReference();
        if (refcount == 0)
        {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
	  cout << "MemoryBlockReference: no more refs, delete MemoryBlock object at "
	       << ((void*)block_) << endl;
#endif

            delete block_;
        }
    }

    void addReference() const
    {
        if (block_) {
#ifdef BZ_DEBUG_LOG_REFERENCES
	  cout << "MemoryBlockReference: reffing block at " << ((void*)block_)
	       << endl;
#endif
            block_->addReference();
        }
        else {
#ifdef BZ_DEBUG_LOG_REFERENCES
            cout << "MemoryBlockReference:: Skipping reference count for data at "<< data_ << endl;
#endif
        }
    };

    int removeReference() const
    {
      if (block_) {
#ifdef BZ_DEBUG_LOG_REFERENCES
	  cout << "MemoryBlockReference: dereffing block at " << ((void*)block_)
	       << endl;
#endif
            return block_->removeReference();
      }
#ifdef BZ_DEBUG_LOG_REFERENCES
        cout << "Skipping reference count for data at "<< data_ << endl;
#endif
        return -1;
    };

    void operator=(const MemoryBlockReference<T_type>&)
    { }
};


BZ_NAMESPACE_END

#include <blitz/memblock.cc>

#endif // BZ_MEMBLOCK_H
