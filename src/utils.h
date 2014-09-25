/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_UTILS
#define H_UTILS

#include "../include/dsrc/Globals.h"

#include <string>

namespace dsrc
{

namespace core
{

template <uint32 _TBitNum>
class TBitMask
{
public:
	static const uint64 Value = ((uint64)1 << (_TBitNum - 1)) | TBitMask<_TBitNum-1>::Value;
};

template <>
class TBitMask<0>
{
public:
	static const uint64 Value = 0;
};


template <uint32 _TNum>
struct TLog2
{
	static const uint32 Value = TLog2< (_TNum >> 1) >::Value + 1;
};

template <>
struct TLog2<1>
{
	static const uint32 Value = 0;
};

template <> struct TLog2<0> {};


template <typename _T>
inline void TSwap(_T& x_, _T& y_)
{
	_T tmp = x_;
	x_ = y_;
	y_ = tmp;
}

template <typename _T>
inline void TFree(_T* p_)
{
	if (p_ != (_T*)0)
		delete p_, p_ = (_T*)0;
}

inline uint32 to_string(uchar* str, uint32 value)
{
	uint32 digits;
	uint32 power = 1;

	if (value == 0)
	{
		str[0] = '0';
		return 1;
	}

	for (digits = 0; digits < 10; ++digits)
	{
		if (value < power)
			break;
		power *= 10;
	}

	power /= 10;
	for (uint32 i = 0; power; ++i, power /= 10)
	{
		int32 d = value / power;
		str[i] = (uchar)('0' + d);
		value -= d * power;
	}

	return digits;
}

inline bool extend_string(uchar *&str, uint32 &size)
{
	uint32 new_size = size * 2;
	uchar *p = new uchar[new_size+1];

	if (!p)
		return false;

	std::copy(str, str+size, p);
	size = new_size;
	delete[] str;
	str = p;

	return true;
}

inline bool extend_string_to(uchar *&str, uint32 &size, uint32 new_size)
{
	if (new_size <= size)
		return true;

	uchar *p = new uchar[new_size+1];

	if (!p)
		return false;

	std::copy(str, str+size, p);

	size = new_size;
	delete[] str;
	str = p;

	return true;
}

inline bool ends_with(const std::string& str_, const std::string& suff_)
{
	return str_.size() >= suff_.size() &&
			str_.compare(str_.size() - suff_.size(), suff_.size(), suff_) == 0;
}

inline uint32 int_log(uint32 x, uint32 base)
{
	uint32 r = 0;

	if (base == 0)
		return 1;
	if (base == 1)
		base++;

	for (uint32 tmp = base; tmp <= x; tmp *= base)
		++r;

	return r;
}

inline uint32 to_num(const uchar *str, uint32 len)
{
	uint32 r = 0;

	for (uint32 i = 0; i < len; ++i)
		r = r * 10 + (str[i] - '0');

	return r;
}

inline bool is_num(const uchar* str_, uint32 len_, uint32& val_)
{
	val_ = 0;
	uint32 i;
	for (i = 0; i < len_; ++i)
	{
		if (str_[i] < '0' || str_[i] > '9')
			break;
		val_ = val_ * 10 + (str_[i] - '0');
	}

	return i == len_ && (len_ == 1 || str_[0] != '0');
}

inline uint32 bit_length(uint64 x)
{
	for (uint32 i = 0; i < 32; ++i)
	{
		if (x < (1ull << i))
			return i;
	}
	return 64;
}

} // namespace core

} // namespace dsrc


#endif
