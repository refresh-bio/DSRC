/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DATASTREAM
#define H_DATASTREAM

#include "../include/dsrc/Globals.h"

namespace dsrc
{

namespace core
{

class IDataStream
{
public:
	virtual ~IDataStream() {}

	virtual int64 PerformIo(uchar* mem_, uint64 size_) = 0;
	virtual void Close() = 0;
};

class IDataStreamReader : public IDataStream
{
public:
	int64 PerformIo(uchar* mem_, uint64 size_)
	{
		return Read(mem_, size_);
	}

	virtual int64 Read(uchar* mem_, uint64 size_) = 0;
};

class IDataStreamWriter : public IDataStream
{
public:
	int64 PerformIo(uchar* mem_, uint64 size_)
	{
		return Write(mem_, size_);
	}

	virtual	int64 Write(const uchar* mem_, uint64 size_) = 0;
};

} // namespace core

} // namespace dsrc

#endif // H_DATASTREAM
