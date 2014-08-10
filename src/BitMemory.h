/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_BITMEMORY
#define H_BITMEMORY

#include "../include/dsrc/Globals.h"

#include <vector>
#include <algorithm>

#include "Buffer.h"
#include "utils.h"


namespace dsrc
{

namespace core
{

class BitMemoryReader
{
public:
	BitMemoryReader(byte *mem_, uint64 size_)
		:	memory(mem_)
		,	size(size_)
		,	position(0)
		,	wordBuffer(0)
		,	wordBufferPos(0)
	{
		ASSERT(mem_ != NULL);
		ASSERT(size_ != 0);
	}

	uint64 Size() const
	{
		return size;
	}

	uint64 Position() const
	{
		return position;
	}

	void SetPosition(uint64 pos_)
	{
		ASSERT(pos_ <= size);
		position = pos_;
	}

	byte* Pointer() const
	{
		return memory;
	}

	// refactor
	//
	uint32 GetBit()
	{
		if (wordBufferPos == 0)
		{
			wordBuffer = GetByte();
			wordBufferPos = 7;
			return (wordBuffer >> 7) & 1;
		}

		return (wordBuffer >> (--wordBufferPos)) & 1;
	}

	uint32 Get2Bits()
	{
		if (wordBufferPos >= 2)
		{
			wordBufferPos -= 2;
			return (wordBuffer >> wordBufferPos) & 3;
		}

		if (wordBufferPos == 0)
		{
			wordBuffer = GetByte();
			wordBufferPos = 6;
			return (wordBuffer >> wordBufferPos) & 3;
		}

		uint32 word = (wordBuffer & 1) << 1;
		wordBuffer = GetByte();
		wordBufferPos = 7;
		word += wordBuffer >> wordBufferPos;
		return word & 3;
	}

	uint32 GetBits(uint32 n_)
	{
		ASSERT(n_ > 0 && n_ < 32);

		uint32 word = 0;
		while (n_)
		{
			if (wordBufferPos == 0)
			{
				wordBuffer = GetByte();
				wordBufferPos = 8;
			}

			if (n_ > wordBufferPos)
			{
				word <<= wordBufferPos;
				word += wordBuffer & BitMask(wordBufferPos);
				n_ -= wordBufferPos;
				wordBufferPos = 0;
			}
			else
			{
				word <<= n_;
				wordBufferPos -= n_;
				word += (wordBuffer >> wordBufferPos) & BitMask(n_);
				break;
			}
		}
		return word;
	}

	byte GetByte()
	{
		ASSERT (position < size);
		return memory[position++];
	}

	void GetBytes(uchar *data, uint32 n_bytes)
	{
		ASSERT(position + n_bytes <= size);

		std::copy(memory + position, memory + position + n_bytes, data);
		position += n_bytes;
	}

	uint32 GetWord()
	{
		uint32 c = GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		return (c << 8) | GetByte();
	}

	uint64 GetDWord()
	{
		uint64 c = GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		return c;
	}

	uint32 Get2Bytes()
	{
		return GetByte() << 8 | GetByte();
	}

	void SkipBytes(uint32 n_)
	{
		ASSERT(position + n_ < size);

		position += n_;
	}

	void FlushInputWordBuffer()
	{
		wordBufferPos = 0;
	}

	void Reset()
	{
		position = 0;
		wordBufferPos = 0;
		wordBuffer = 0;
	}

private:
	static const uint32 WordBufferSize = 8;

	byte* memory;
	uint64 size;
	uint64 position;

	uint32 wordBuffer;
	uint32 wordBufferPos;

	BitMemoryReader(const BitMemoryReader& )
	{}

	BitMemoryReader& operator= (const BitMemoryReader& )
	{ return *this; }


	uint32 BitMask(uint32 n_)
	{
		ASSERT(n_ < 32);
		//return ((uint32)0xFFFFFFFF >> (32 - n_));
		return ((uint32)1 << n_) - 1;
	}
};


class BitMemoryWriter
{
public:
	static const uint32 DefaultBufferSize = 1 << 20;

	BitMemoryWriter(uint32 bufferSize_ = DefaultBufferSize)
		:	buffer(NULL)
		,	position(0)
		,	wordBuffer(0)
		,	wordBufferPos(0)
		,	hasOwnership(true)
	{
		buffer = new Buffer(bufferSize_);
		memory = buffer->Pointer();
		size = buffer->Size();
	}

	BitMemoryWriter(Buffer& buffer_)
		:	buffer(NULL)
		,	position(0)
		,	wordBuffer(0)
		,	wordBufferPos(0)
		,	hasOwnership(false)
	{
		buffer = &buffer_;
		memory = buffer->Pointer();
		size = buffer->Size();
	}

	~BitMemoryWriter()
	{
		ASSERT(buffer != NULL);
		if (hasOwnership)
		{
			delete buffer;
		}
	}

	uint64 Size() const
	{
		return size;
	}

	uint64 Position() const
	{
		return position;
	}

	void SetPosition(uint64 pos_)
	{
		ASSERT(pos_ <= size);
		position = pos_;
	}

	byte* Pointer() const
	{
		return memory;
	}

	template <typename _T>
	void PutBit(_T b_)
	{
		if (wordBufferPos < WordBufferSize)
		{
			wordBuffer <<= 1;
			wordBuffer += ((uint32)b_) & 1;
			++wordBufferPos;
		}
		else
		{
			PutWord(wordBuffer);
			wordBufferPos = 1;
			wordBuffer = ((uint32)b_) & 1;
		}
	}

	void Put2Bits(uint32 word_)
	{
		ASSERT(word_ <= 3);
		word_ &= 3;

		if (wordBufferPos + 2 <= WordBufferSize)
		{
			wordBuffer <<= 2;
			wordBuffer += word_;
			wordBufferPos += 2;
		}
		else if (wordBufferPos == WordBufferSize)
		{
			PutWord(wordBuffer);
			wordBufferPos = 2;
			wordBuffer = word_;
		}
		else
		{
			wordBuffer <<= 1;
			wordBuffer += word_ >> 1;
			PutWord(wordBuffer);
			wordBuffer = word_ & 1;
			wordBufferPos = 1;
		}
	}

	void PutBits(uint32 word_, uint32 n_)
	{
		//ASSERT(word <= n_bit_mask[n_bits]);
		ASSERT(n_ > 0);
		word_ &= BitMask(n_);

		uint32 rest = WordBufferSize - wordBufferPos;
		if (n_ >= rest)
		{
			n_ -= rest;
			wordBuffer <<= rest;
			wordBuffer += word_ >> n_;
			wordBufferPos = 0;
			PutWord(wordBuffer);
			wordBuffer = 0;
		}

		wordBuffer <<= n_;
		wordBuffer += word_ & BitMask(n_);
		wordBufferPos += n_;
	}

	void PutByte(byte b_)
	{
		if (position >= size)
		{
			ExtendBuffer(size + (size >> 2));
		}

		memory[position++] = b_;
	}

	void Put2Bytes(uint32 word_)
	{
		PutByte((uchar)(word_ >> 8));
		PutByte((uchar)(word_ & 0xFF));
	}

	void PutBytes(const byte *data_, uint32 n_)
	{
		if (position + n_ > size)
		{
			ExtendBuffer((uint32)((position + n_) * 1.5f));
		}
		std::copy(data_, data_ + n_, memory + position);
		position += n_;
	}

	void PutDWord(uint64 data_)
	{
		PutByte(data_ >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
	}

	void PutWord(uint32 data_)
	{
		PutByte(data_ >> 24);
		PutByte((data_ >> 16) & 0xFF);
		PutByte((data_ >> 8) & 0xFF);
		PutByte(data_ & 0xFF);
	}

	void FlushFullWordBuffer()
	{
		PutWord(wordBuffer);

		wordBuffer    = 0;
		wordBufferPos = 0;
	}

	void FlushPartialWordBuffer()
	{
		wordBuffer <<= (32 - wordBufferPos) & 7;

		if (wordBufferPos > 24)
			PutByte(wordBuffer >> 24);
		if (wordBufferPos > 16)
			PutByte((wordBuffer >> 16) & 0xFF);
		if (wordBufferPos > 8)
			PutByte((wordBuffer >> 8) & 0xFF);
		if (wordBufferPos > 0)
			PutByte(wordBuffer & 0xFF);

		wordBuffer = 0;
		wordBufferPos = 0;
	}

	void Flush()
	{
		FlushPartialWordBuffer();
	}

	void SwapFlush(Buffer &buf_, uint64& size_)
	{
		ASSERT(!hasOwnership);
		ASSERT(buffer != 0 && position != 0);

		FlushPartialWordBuffer();
		size_ = position;

		buffer->Swap(buf_);
		memory = buffer->Pointer();
		size = buffer->Size();
		position = 0;
	}

	void Reset()
	{
		position = 0;
		wordBufferPos = 0;
		wordBuffer = 0;
	}


private:
	static const uint32 WordBufferSize = 32;

	Buffer*	buffer;
	byte*	memory;
	uint64	size;
	uint64	position;

	uint32 wordBuffer;
	uint32 wordBufferPos;

	bool	hasOwnership;

	BitMemoryWriter(const BitMemoryWriter& )
		:	buffer(0)
	{}

	BitMemoryWriter& operator= (const BitMemoryWriter& )
	{ return *this; }


	uint32 BitMask(uint32 n_)
	{
		ASSERT(n_ < 32);

		return ((uint32)1 << n_) - 1;
	}

	void ExtendBuffer(uint32 newSize_)
	{
		ASSERT(newSize_ > 0);

		buffer->Extend(newSize_, true);
		ASSERT(newSize_ == buffer->Size());

		memory = buffer->Pointer();
		size = buffer->Size();
	}
};

} // namespace core

} // namespace dsrc


#endif
