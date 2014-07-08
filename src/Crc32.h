/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#ifndef CRC32_H
#define CRC32_H

#include "../include/dsrc/Globals.h"

#include <vector>
#include <string>

namespace dsrc
{

namespace core
{

class Crc32Hasher	// CRC32 LSB
{
private:
	uint32 polynomial;
	uint32 crc;
	uint32 lookup_table[256];

	void FillLookupTable()
	{
		lookup_table[0] = 0;
		for (uint32 i = 1; i < 256; ++i)
		{
			uint32 h = i;
			for (uint32 j = 0; j < 8; ++j)
			{
				if (h & 1)
				{
					h = polynomial ^ (h >> 1);
				}
				else
				{
					h >>= 1;
				}
			}
			lookup_table[i] = h;
		}
	}

public:
	Crc32Hasher(uint32 polynomial_ = 0xEDB88320, uint32 seed_ = 0xFFFFFFFF)
		:	polynomial(polynomial_)
		,	crc(seed_)
	{
		FillLookupTable();
	}

	void UpdateCrc(uchar c_)
	{
		crc = (crc >> 8) ^ lookup_table[(c_ ^ crc) & 0xFF];
	}

	void UpdateCrc(const uchar* str_, uint32 len_)
	{
		for (uint32 i = 0; i < len_; ++i)
		{
			UpdateCrc(str_[i]);
		}
	}

	uint32 GetHash() const
	{
		return crc ^ 0xFFFFFFFF;
	}

	void Reset(uint32 polynomial_ = 0xEDB88320, uint32 seed_ = 0xFFFFFFFF)
	{
		if (polynomial != polynomial_)
		{
			polynomial = polynomial_;
			FillLookupTable();
		}
		crc = seed_;
	}

	uint32 ComputeHash(const uchar* arr_, uint32 len_)
	{
		ASSERT(arr_ != NULL);

		Reset();
		for (uint32 i = 0; i < len_; ++i)
		{
			UpdateCrc(arr_[i]);
		}
		return GetHash();
	}

	uint32 ComputeHash(const std::string& str_)
	{
		return Crc32Hasher::ComputeHash((const uchar*)str_.c_str(), (uint32)str_.length());
	}
};

} // namespace core

} // namespace dsrc

#endif
