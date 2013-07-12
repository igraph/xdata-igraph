/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef BLOCKMODEL_IO_HPP
#define BLOCKMODEL_IO_HPP

#include <iostream>
#include <string>
#include <block/blockmodel.h>

/// Supported formats in the IO subsystem
typedef enum {
    FORMAT_NULL, FORMAT_PLAIN, FORMAT_JSON
} Format;

/***************************************************************************/

/// Abstract templatized reader class
/**
 * Specializations of this class must implement reading (deserializing)
 * some kind of an object in some format from an input stream.
 */
template <typename T>
class Reader {
protected:
    /// The name of the file the reader read from the last time (if known)
    std::string m_filename;

public:
    /// Returns the name of the file the reader read from the last time
    const std::string getFilename() const {
        return m_filename;
    }

    /// Reads into the given object from the given stream
    virtual void read(T* model, std::istream& is) = 0;
};

/// Reader for objects in plain text format
template <typename T>
class PlainTextReader : public Reader<T> {
public:
    /// Writes the given object to the given stream
    virtual void read(T* model, std::istream& is);
};

/// Reader for blockmodels in plain text format
template <>
class PlainTextReader<Blockmodel> : public Reader<Blockmodel> {
public:
    /// Initializes the given blockmodel from the given stream
    virtual void read(Blockmodel* model, std::istream& is);
};

/***************************************************************************/

/// Abstract templatized writer class
/**
 * Specializations of this class must implement writing (serializing)
 * some kind of an object in some format to an output stream.
 */
template <typename T>
class Writer {
public:
    /// Writes the given object to the given stream
    virtual void write(const T& model, std::ostream& os) {
        write(&model, os);
    }

    /// Writes the given object to the given stream
    virtual void write(const T* model, std::ostream& os) = 0;
};

/// Null writer that does nothing
template <typename T>
class NullWriter : public Writer<T> {
public:
    /// Does nothing.
    virtual void write(const T*, std::ostream&) {}
};

/// Writer for objects in plain text format
template <typename T>
class PlainTextWriter : public Writer<T> {
public:
    /// Writes the given object to the given stream
    virtual void write(const T* model, std::ostream& os);
};

/// Writer for objects in JSON format
template <typename T>
class JSONWriter : public Writer<T> {
public:
    /// Writes the given object to the given stream
    virtual void write(const T* model, std::ostream& os);
};

#endif

