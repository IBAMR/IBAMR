/**
 * @file MyHeaps.h
 * @author Andrea Tagliasacchi
 * @data 26 March 2008
 * @copyright (c) Andrea Tagliasacchi - All rights reserved
 */

#ifndef MYHEAPS_H_
#define MYHEAPS_H_

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <exception> // general exception
#include <iostream>
#include <stdexcept> // out_of_range
#include <vector>

// macros for navigation in the hard coded binary tree
#define PARENT(pos) ((pos - 1) >> 1) // equivalent to floor(pos/2)
#define LEFT(pos) ((pos << 1) + 1)   // equivalent to pos*2 + 1
#define RIGHT(pos) ((pos << 1) + 2)  // equivalent to pos*2 + 2

/// EXCEPTION
class HeapEmptyException : public std::out_of_range
{
public:
    HeapEmptyException(const std::string& message) : std::out_of_range(message)
    {
        ;
    }
};
class InvalidKeyIncreaseException : public std::out_of_range
{
public:
    InvalidKeyIncreaseException(const std::string& message) : std::out_of_range(message)
    {
        ;
    }
};
class InvalidIndexException : public std::out_of_range
{
public:
    InvalidIndexException(const std::string& message) : std::out_of_range(message)
    {
        ;
    }
};

/**
 * This class provides a back-inxedex heap structure where indexes of
 * elements already in the heap are kept updated to allow for random access
 * update of elements (done automatically in push if element with
 * "idx" is already contained in the heap )
 *
 * Refer to the following textbook for details:
 * @book{cormen1990ia,
 *   title={{Introduction to algorithms}},
 *   author={Cormen, T.T. and Leiserson, C.E. and Rivest, R.L.},
 *   year={1990},
 *   publisher={MIT Press Cambridge, MA, USA}
 * }
 */
template <class Tkey>
class MaxHeap
{
private:
    /// root is assumed to be at end of the vector
    std::vector<std::pair<Tkey, int> > heap;
    /**
     *  maintain a list of back indexes.
     *  * -1           not in heap
     *  * other 	   index that point to cell in vector heap
     */
    std::vector<int> backIdx;
    /**
     * If useBackIdx==false it means that the current structure
     * is not making use of a backindexed heap. Thus, no update
     * is available
     */
    bool useBackIdx;

public:
    /// Simple constructor with NO cross updates
    MaxHeap()
    {
        useBackIdx = false;
    }
    /// back indexes constructor used for cross updates
    MaxHeap(int Nindex)
    {
        // initialize the back indexes with pseudo-null pointers
        backIdx.resize(Nindex, -1);
        useBackIdx = true;
    }

    /// pushes a new value in the heap
    void push(Tkey key, int index)
    {
        // cout << "pushing " << index << endl;
        if (useBackIdx && index >= (int)backIdx.size())
            throw InvalidIndexException(
                "the index in the push must be smaller than the maximal allowed index (specified in constructor)");

        // If key is not in backindexes or there is no backindexes AT ALL.... complete push (no update)
        if (!useBackIdx)
        {
            // add to the back of the vector
            heap.emplace_back(std::make_pair(key, index));
            // recursive call to increase key
            heapIncreaseKey(static_cast<int>(heap.size()) - 1, key);
        }
        else
        {
            if (backIdx[index] == -1)
            {
                // add to the back of the vector
                heap.emplace_back(std::make_pair(key, index));
                // initially point to back
                backIdx[index] = static_cast<int>(heap.size()) - 1;
                // recursive call to increase key
                heapIncreaseKey(static_cast<int>(heap.size()) - 1, key);
                // USE STL STUFF
                // push_heap(heap.begin(),heap.end());
            }
            // update push (a key exists)
            else
            {
                heapIncreaseKey(backIdx[index], key);
            }
        }
    }

    /// return a constant reference to the MINIMAL KEY element stored in the head of the heap
    const std::pair<Tkey, int>& top()
    {
        if (heap.empty())
        {
            std::cout << "Impossible to get top element, empty heap";
            std::cin.get();
            std::exit(0);
        }
        else
            return heap[0];
    }

    /// removes the top element of the queue (minimal)
    void pop()
    {
        if (heap.size() < 1)
        { // a.k.a. heap.empty()

            std::cout << "Heap underflow";
            std::cin.get();
            std::exit(0);
        }
        // overwrite top with tail element
        heap[0] = heap.back();

        // USE STL FUNCTIONALITIES (NOT ALLOW BACKINDEXs)
        // pop_heap(heap.begin(), heap.end());

        // shorten the vector
        heap.pop_back();

        // start heapify from root
        maxHeapify(0);
    }

    /// returns the size of the heap
    int size()
    {
        return static_cast<int>(heap.size());
    }

    /// check for emptyness
    bool empty()
    {
        return heap.empty();
    }

    /// check recursively if the substructures is correct using STL provided algorithm
    bool verifyHeap()
    {
        return std::is_heap(heap.begin(), heap.end());
    }

private:
    /// check and applies MaxHeap Correctness down the subtree with index "currIdx"
    void maxHeapify(int currIdx)
    {
        unsigned int leftIdx = LEFT(currIdx);
        unsigned int rightIdx = RIGHT(currIdx);

        // decide if and where ta swap, left or right, then swap
        // current is the best choice (defalut)
        int largestIdx;

        // is left a better choice? (exists an invalid placed bigger value on the left side)
        if (leftIdx < heap.size() && heap[leftIdx].first > heap[currIdx].first)
            largestIdx = leftIdx;
        else
            largestIdx = currIdx;

        // is right a better choice? (exists an invalid placed bigger value on the right side)
        if (rightIdx < heap.size() && heap[rightIdx].first > heap[largestIdx].first) largestIdx = rightIdx;

        // a better choice exists?
        if (largestIdx != currIdx)
        {
            // swap elements
            swap(currIdx, largestIdx);

            // recursively call this function on alterated subtree
            maxHeapify(largestIdx);
        }
    }

    /// swap the content of two elements in position pos1 and pos2
    void swap(int pos1, int pos2)
    {
        assert(!heap.empty());
        assert(pos1 >= 0 && pos1 < (int)heap.size());
        assert(pos2 >= 0 && pos2 < (int)heap.size());

        // update backindexes
        if (useBackIdx)
        {
            backIdx[heap[pos1].second] = pos2;
            backIdx[heap[pos2].second] = pos1;
        }

        // update heap
        std::pair<Tkey, int> temp = heap[pos1];
        heap[pos1] = heap[pos2];
        heap[pos2] = temp;
    }

    /// propagates the correctness (in heap sense) down from a vertex currIdx
    void heapIncreaseKey(int currIdx, Tkey key)
    {
        // check if given key update is actually an increase
        if (key < heap[currIdx].first)
            throw InvalidKeyIncreaseException("In MaxHeaps only increase key updates are legal");

        // update value with current key
        heap[currIdx].first = key;

        // traverse the tree up making necessary swaps
        int parentIdx = PARENT(currIdx);
        while (currIdx > 0)
        {
            if (heap[parentIdx].first < heap[currIdx].first)
            {
                // make swap
                swap(currIdx, parentIdx);
                // move up
                currIdx = parentIdx;
                parentIdx = PARENT(currIdx);
            }
            else
            {
                break;
            }
        }
    }

    /// print an internal representation of the heap (debug purposes)
    void print()
    {
        std::cout << "idxs";
        for (int i = 0; i < size(); i++) std::cout << " " << heap[i].second << " ";
        std::cout << std::endl;

        std::cout << "csts";
        for (int i = 0; i < size(); i++) std::cout << " " << heap[i].first << " ";
        std::cout << std::endl;

        //		cout << "";
        //		for ( int i=0; i < size(); i++)
        //			cout << heap[i].first << " in off: " << backIdx[heap[i].first] << ", ";
        //		cout << endl;

        std::cout << std::endl;
    }
};

/**
 * This class provides a back-inxedex heap (MinHeap) structure where indexes of
 * elements already in the heap are kept updated to allow for random access
 * update of elements (done automatically in push if element with
 * "idx" is already contained in the heap )
 *
 * Refer to the following textbook for details:
 * @book{cormen1990ia,
 *   title={{Introduction to algorithms}},
 *   author={Cormen, T.T. and Leiserson, C.E. and Rivest, R.L.},
 *   year={1990},
 *   publisher={MIT Press Cambridge, MA, USA}
 * }
 */
template <class Tkey>
class MinHeap
{
private:
    /// root is assumed to be at end of the vector
    std::vector<std::pair<Tkey, int> > heap;
    /**
     *  maintain a list of back indexes.
     *  * -1           not in heap
     *  * other 	   index that point to cell in vector heap
     */
    std::vector<int> backIdx;
    /**
     * If useBackIdx==false it means that the current structure
     * is not making use of a backindexed heap. Thus, no update
     * is available
     */
    bool useBackIdx;

public:
    /// back indexes constructor used for cross updates
    MinHeap(int Nindex)
    {
        // initialize the back indexes with pseudo-null pointers
        backIdx.resize(Nindex, -1);
        useBackIdx = true;
    }
    /// Simple constructor with NO cross updates
    MinHeap()
    {
        useBackIdx = false;
    }

    /// pushes a new value in the heap
    void push(Tkey key, int index)
    {
        // cout << "pushing " << index << endl;
        if (useBackIdx && index >= (int)backIdx.size())
            throw InvalidIndexException(
                "the index in the push must be smaller than the maximal allowed index (specified in constructor)");

        // If key is not in backindexes or there is no backindexes AT ALL.... complete push (no update)
        if (!useBackIdx)
        {
            // add to the back of the vector
            heap.push_back(std::make_pair(key, index));
            // recursive call to increase key
            heapDecreaseKey(static_cast<int>(heap.size()) - 1, key);
        }
        else
        {
            if (useBackIdx || backIdx[index] == -1)
            {
                // add to the back of the vector
                heap.push_back(std::make_pair(key, index));
                // initially point to back
                backIdx[index] = static_cast<int>(heap.size()) - 1;
                // recursive call to increase key
                heapDecreaseKey(static_cast<int>(heap.size()) - 1, key);
                // USE STL STUFF
                // push_heap(heap.begin(),heap.end());
            }
            // update push (a key exists)
            else
            {
                heapDecreaseKey(backIdx[index], key);
            }
        }
    }

    /// return a constant reference to the MINIMAL KEY element stored in the head of the heap
    const std::pair<Tkey, int>& top()
    {
        if (heap.empty())
        {
            std::cout << "Impossible to get top element, empty heap";
            std::cin.get();
            std::exit(0);
        }
        else
            return heap[0];
    }

    /// removes the top element of the queue (minimal)
    void pop()
    {
        if (heap.size() < 1) // a.k.a. heap.empty()
        {
            std::cout << "heap underflow";
            std::cin.get();
        }
        // overwrite top with tail element
        heap[0] = heap.back();

        // USE STL FUNCTIONALITIES (NOT ALLOW BACKINDEXs)
        // pop_heap(heap.begin(), heap.end());

        // shorten the vector
        heap.pop_back();

        // start heapify from root
        minHeapify(0);
    }

    /// returns the size of the heap
    int size()
    {
        return heap.size();
    }

    /// check for emptyness
    bool empty()
    {
        return heap.empty();
    }

    // this does not work, how do you provide a new ordering function to is_heap??
    /// check recursively if the substructures is correct using STL provided algorithm
    // bool verifyHeap( ){
    //	return std::__is_heap(heap.begin(), heap.end() );
    //}

    /// computes full heap sort and returns the corresponding indexing structure
    /// Requires the indexes to be allocated already.
    void heapsort(std::vector<int>& indexes)
    {
        // until empty... keep popping
        int i = 0;
        while (empty() == false)
        {
            std::pair<Tkey, int> t = top();
            pop();
            indexes[i++] = t.second;
        }
    }

private:
    /// check and applies MaxHeap Correctness down the subtree with index "currIdx"
    void minHeapify(int currIdx)
    {
        unsigned int leftIdx = LEFT(currIdx);
        unsigned int rightIdx = RIGHT(currIdx);

        // decide if and where ta swap, left or right, then swap
        // current is the best choice (defalut)
        int smallerIdx;

        // is left a better choice? (exists an invalid placed smaller value on the left side)
        if (leftIdx < heap.size() && heap[leftIdx].first < heap[currIdx].first)
            smallerIdx = leftIdx;
        else
            smallerIdx = currIdx;

        // is right a better choice? (exists an invalid placed smaller value on the right side)
        if (rightIdx < heap.size() && heap[rightIdx].first < heap[smallerIdx].first) smallerIdx = rightIdx;

        // a better choice exists?
        if (smallerIdx != currIdx)
        {
            // swap elements
            swap(currIdx, smallerIdx);

            // recursively call this function on alterated subtree
            minHeapify(smallerIdx);
        }
    }

    /// swap the content of two elements in position pos1 and pos2
    void swap(int pos1, int pos2)
    {
        assert(!heap.empty());
        assert(pos1 >= 0 && pos1 < (int)heap.size());
        assert(pos2 >= 0 && pos2 < (int)heap.size());

        // update backindexes
        if (useBackIdx)
        {
            backIdx[heap[pos1].second] = pos2;
            backIdx[heap[pos2].second] = pos1;
        }

        // update heap
        std::pair<Tkey, int> temp = heap[pos1];
        heap[pos1] = heap[pos2];
        heap[pos2] = temp;
    }

    /// propagates the correctness (in heap sense) down from a vertex currIdx
    void heapDecreaseKey(int currIdx, Tkey key)
    {
        // check if given key update is actually an increase
        if (key > heap[currIdx].first)
            throw InvalidKeyIncreaseException("In MinHeaps only decrease in key updates are legal");

        // update value with current key
        heap[currIdx].first = key;

        // traverse the tree up making necessary swaps
        int parentIdx = PARENT(currIdx);
        while (currIdx > 0)
        {
            if (heap[parentIdx].first > heap[currIdx].first)
            {
                // make swap
                swap(currIdx, parentIdx);
                // move up
                currIdx = parentIdx;
                parentIdx = PARENT(currIdx);
            }
            else
            {
                break;
            }
        }
    }

    /// print an internal representation of the heap (debug purposes)
public:
    void print()
    {
        std::cout << "idxs";
        for (int i = 0; i < size(); i++) std::cout << " " << heap[i].second << " ";
        std::cout << std::endl;

        std::cout << "csts";
        for (int i = 0; i < size(); i++) std::cout << " " << heap[i].first << " ";
        std::cout << std::endl;

        //		cout << "";
        //		for ( int i=0; i < size(); i++)
        //			cout << heap[i].first << " in off: " << backIdx[heap[i].first] << ", ";
        //		cout << endl;

        std::cout << std::endl;
    }
};

#endif /*MYHEAPS_H_*/
