/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef BUFFER_H
#define BUFFER_H

#include "../include/dsrc/Globals.h"

#include <algorithm>

#include "utils.h"


namespace dsrc
{

namespace core
{


class Buffer
{
public:
	Buffer(uint64 size_)
		:	size(size_)
		,	ownsMemory(true)
	{
		ASSERT(size_ != 0);

		buffer = new byte[size_];
	}

	Buffer(byte* mem_, uint64 size_)
		:	buffer(mem_)
		,	size(size_)
		,	ownsMemory(false)
	{
		ASSERT(size_ != 0);
		ASSERT(mem_ != NULL);
	}

	~Buffer()
	{
		if (ownsMemory)
			delete buffer;
	}

	uint64 Size() const
	{
		return size;
	}

	byte* Pointer() const
	{
		return (byte*)buffer;
	}

	void Extend(uint64 size_, bool copy_ = false)
	{
		ASSERT(ownsMemory);		// TODO: should be an error thrown here

		if (size > size_)
			return;

		byte* p = new byte[size_];

		if (copy_)
			std::copy(buffer, buffer + size, p);

		delete[] buffer;

		buffer = p;
		size = size_;
	}

	void Swap(Buffer& b)
	{
		TSwap(b.buffer, buffer);
		TSwap(b.size, size);
		TSwap(b.ownsMemory, ownsMemory);
	}

private:
	Buffer(const Buffer& ) {}
	Buffer& operator= (const Buffer& )
	{ return *this; }

	byte* buffer;
	uint64 size;
	bool ownsMemory;
};


struct DataChunk
{
	static const uint64 DefaultBufferSize = 1 << 20;		// 1 << 22

	Buffer data;
	uint64 size;

	DataChunk(uint64 bufferSize_ = DefaultBufferSize)
		:	data(bufferSize_)
		,	size(0)
	{}

	DataChunk(byte* buffer_, uint64 bufferSize_)
		:	data(buffer_, bufferSize_)
		,	size(0)
	{}

	void Reset()
	{
		size = 0;
	}
};

} // namespace core

} // namespace dsrc


#endif // BUFFER_H
