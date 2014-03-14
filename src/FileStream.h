/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
#ifndef H_FILESTREAM
#define H_FILESTREAM

#include "../include/dsrc/Globals.h"

#include "DataStream.h"

namespace dsrc
{

namespace core
{

class IFileStream
{
public:
	IFileStream();
	virtual ~IFileStream();

protected:
	struct FileStreamImpl;
	FileStreamImpl* impl;
};

class FileStreamReader : public IDataStreamReader, public IFileStream
{
public:
	FileStreamReader(const std::string& fileName_);
	~FileStreamReader();

	void Close();

	virtual int64 Read(uchar* mem_, uint64 size_);
};

class FileStreamReaderExt : public FileStreamReader
{
public:
	FileStreamReaderExt(const std::string& fileName_);

	void SetPosition(uint64 pos_);

	uint64 Position() const
	{
		return position;
	}

	uint64 Size() const
	{
		return size;
	}

	virtual int64 Read(uchar* mem_, uint64 size_);

private:
	uint64 size;
	uint64 position;
};


class FileStreamWriter : public IDataStreamWriter, public IFileStream
{
public:
	FileStreamWriter(const std::string& fileName_);
	~FileStreamWriter();

	void Close();

	int64 PerformIo(uchar* mem_, uint64 size_)
	{
		return Write(mem_, size_);
	}

	virtual	int64 Write(const uchar* mem_, uint64 size_);
};

class FileStreamWriterExt : public FileStreamWriter
{
public:
	FileStreamWriterExt(const std::string& fileName_);

	void SetPosition(uint64 pos_);

	uint64 Position() const
	{
		return position;
	}

	virtual	int64 Write(const uchar* mem_, uint64 size_);

private:
	uint64 position;
};

} // namespace core

} // namespace dsrc

#endif // H_FILESTREAM
