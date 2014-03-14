/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_STDSTREAM
#define H_STDSTREAM

#include "../include/dsrc/Globals.h"

#include "DataStream.h"

namespace dsrc
{

namespace core
{

class StdStreamReader : public IDataStreamReader
{
public:
	StdStreamReader()
	{}

	~StdStreamReader()
	{}

	void Close()
	{}

	int64 Read(uchar* mem_, uint64 size_);
};

class StdStreamWriter : public IDataStreamWriter
{
public:
	StdStreamWriter()
	{}

	~StdStreamWriter()
	{}

	void Close()
	{}

	int64 Write(const uchar* mem_, uint64 size_);
};

} // namespace core

} // namespace dsrc

#endif // H_STDSTREAM
