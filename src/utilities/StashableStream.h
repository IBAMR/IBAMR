#ifndef included_StashableStream
#define included_StashableStream

// Filename: StashableStream.h
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <07.Apr.2007 19:27:43 griffith@box221.cims.nyu.edu>

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * Class StashableStream implements a message buffer of fixed size
 * used by the communication routines.  It implements the
 * SAMRAI::tbox::AbstractStream interface.  Class StashableStream can
 * packs and unpacks message streams via straight-forward byte
 * copying.
 *
 * \note This class will not work on heterogeneous machines, which
 * require a machine-independent mechanism such as XDR.
 */
class StashableStream
    : public SAMRAI::tbox::AbstractStream
{
public:
    enum StreamMode { Read, Write };

    /*!
     * Create a message stream of the specified size in bytes and the
     * stream mode (one of StashableStream::Read or
     * StashableStream::Write).
     */
    StashableStream(
        const int bytes,
        const StreamMode mode);

    /*!
     * Create a message stream with the specified buffer and stream
     * mode (one of StashableStream::Read or StashableStream::Write).
     */
    StashableStream(
        const void* const buffer,
        const int bytes,
        const StreamMode mode);

    /*!
     * Virtual destructor for a message stream.
     */
    virtual ~StashableStream();

    /*!
     * Return a pointer to the start of the message buffer.
     */
    void* getBufferStart();

    /*!
     * Return the current size of the buffer in bytes.
     */
    int getCurrentSize() const;

    /*!
     * Return the current index into the buffer.
     */
    int getCurrentIndex() const;

    /*!
     * Set the current index into the buffer.  Further
     * packing/unpacking will begin at this new location.
     */
    void setCurrentIndex(
        const int index);

    /*!
     * Reset the index to the beginning of the buffer.  This is the
     * same as setting the buffer index to zero via setCurrentIndex().
     */
    void resetIndex();

    /*!
     * @name Boolean Stream Primitives
     * Pack and unpack booleans into and out of the message stream.
     */
    //@{
    /// Pack a single bool into the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator<<(const bool& data);
    /// Remove a single bool from the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator>>(bool& data);
    /// Pack an array of bools into the message stream.
    virtual void pack(const bool* data, const int n=1);
    /// Remove an array of bools from the message stream.
    virtual void unpack(bool* data, const int n=1);
    //@}

    /*!
     * @name Character Stream Primitives
     * Pack and unpack chars into and out of the message stream.
     */
    //@{
    /// Pack a single char into the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator<<(const char& data);
    /// Remove a single char from the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator>>(char& data);
    /// Pack an array of chars into the message stream.
    virtual void pack(const char* data, const int n=1);
    /// Remove an array of chars from the message stream.
    virtual void unpack(char* data, const int n=1);
    //@}

    /*!
     * @name Double Complex Stream Primitives
     * Pack and unpack double complex into and out of the message stream.
     */
    //@{
    /// Pack a single double complex into the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator<<(const dcomplex& data);
    /// Remove a single double complex from the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator>>(dcomplex& data);
    /// Pack an array of double complex into the message stream.
    virtual void pack(const dcomplex* data, const int n=1);
    /// Remove an array of double complex from the message stream.
    virtual void unpack(dcomplex* data, const int n=1);
    //@}

    /*!
     * @name Double Stream Primitives
     * Pack and unpack doubles into and out of the message stream.
     */
    //@{
    /// Pack a single double into the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator<<(const double& data);
    /// Remove a single double from the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator>>(double& data);
    /// Pack an array of doubles into the message stream.
    virtual void pack(const double* data, const int n=1);
    /// Remove an array of doubles from the message stream.
    virtual void unpack(double* data, const int n=1);
    //@}

    /*!
     * @name Float Stream Primitives
     * Pack and unpack floats into and out of the message stream.
     */
    //@{
    /// Pack a single float into the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator<<(const float& data);
    /// Remove a single float from the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator>>(float& data);
    /// Pack an array of floats into the message stream.
    virtual void pack(const float* data, const int n=1);
    /// Remove an array of floats from the message stream.
    virtual void unpack(float* data, const int n=1);
    //@}

    /*!
     * @name Integer Stream Primitives
     * Pack and unpack integers into and out of the message stream.
     */
    //@{
    /// Pack a single integer into the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator<<(const int& data);
    /// Remove a single integer from the message stream.
    virtual SAMRAI::tbox::AbstractStream& operator>>(int& data);
    /// Pack an array of integers into the message stream.
    virtual void pack(const int* data, const int n=1);
    /// Remove an array of integers from the message stream.
    virtual void unpack(int* data, const int n=1);
    //@}

    /*!
     * Print out internal class data for debugging.
     */
    virtual void printClassData(
        ostream& os) const;

private:
    /*!
     * @brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    StashableStream();

    /*!
     * @brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    StashableStream(
        const StashableStream& from);

    /*!
     * @brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    StashableStream& operator=(
        const StashableStream& that);

    /*!
     * \brief Return a pointer to buffer space and advance internal
     * pointers to reflect the allocated buffers space.
     */
    void* getPointerAndAdvanceCursor(
        const int bytes);

    /*!
     * \brief Pack the specified data to the buffer.
     */
    template<typename T>
    void __pack(
        const T* const m_data,
        unsigned m_bytes);

    /*!
     * \brief Unpack the specified data to the buffer.
     */
    template<typename T>
    void __unpack(
        T* const m_data,
        unsigned m_bytes);

    /*
     * The size of the buffer.
     */
    const int d_buffer_size;

    /*
     * The current size of the buffer, i.e., the number of bytes in
     * the buffer which are currently in use.
     */
    int d_current_size;

    /*
     * The index of the first free element in the buffer.
     */
    int d_buffer_index;

    /*
     * The buffer.
     */
    char* const d_buffer;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/StashableStream.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StashableStream
