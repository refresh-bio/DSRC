/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_SYMBOLCODERRC
#define H_SYMBOLCODERRC

#include "../include/dsrc/Globals.h"

#include "RangeCoder.h"

namespace dsrc
{

namespace comp
{

template <uint32 _TMaxSymbolCount>
class TSymbolCoderRC
{
public:
	typedef uint16 StatType;
	static const uint32 MaxSymbolCount = _TMaxSymbolCount > 0 ? _TMaxSymbolCount : 1;

	TSymbolCoderRC()
	{
		std::fill(stats, stats + MaxSymbolCount, 1);
	}

	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		ASSERT(sym_ < MaxSymbolCount);

		uint32 acc = Accumulate();
		uint32 loEnd = 0;

		for (uint32 i = 0; i < sym_; ++i)
			loEnd += stats[i];

		rc_.EncodeFrequency(stats[sym_], loEnd, acc);

		stats[sym_] += StepSize;
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 acc = Accumulate();
		uint32 cul = rc_.GetCumulativeFreq(acc);

		uint32 idx, hiEnd;
		for (idx = 0, hiEnd = 0; (hiEnd += stats[idx]) <= cul; ++idx)
			;
		hiEnd -= stats[idx];

		rc_.UpdateFrequency(stats[idx], hiEnd, acc);
		stats[idx] += StepSize;
		return idx;
	}

private:
	static const StatType StepSize = 2;
	static const uint32 MaxAccumulatedValue = (1<<16) - MaxSymbolCount*StepSize;

	void Rescale()
	{
		for (uint32 i = 0; i < MaxSymbolCount; ++i)
			stats[i] -= stats[i] >> 1;		// no '>>=' to avoid reducing stats to 0
	}

	uint32 Accumulate()
	{
		uint32 acc = 0;
		for (uint32 i = 0; i < MaxSymbolCount; ++i)
			acc += stats[i];

		if (acc >= MaxAccumulatedValue)
		{
			acc = 0;
			Rescale();
			for (uint32 i = 0; i < MaxSymbolCount; ++i)
				acc += stats[i];
		}

		return acc;
	}

	StatType stats[MaxSymbolCount];			// can be assumed to be uint16
};

class SymbolCoderRC
{
public:
	typedef uint16 StatType;
	static const uint32 MaxSymbolCount = 256;

	SymbolCoderRC()
		:	symbolCount(MaxSymbolCount)
	{
		std::fill(stats, stats + MaxSymbolCount, 1);
	}

	void SetSymbolCount(uint32 count_)
	{
		ASSERT(count_ > 0);
		ASSERT(count_ <= MaxSymbolCount);

		symbolCount = count_;
	}

	void Clear()
	{
		std::fill(stats, stats + MaxSymbolCount, 1);
	}

	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		ASSERT(sym_ < symbolCount);

		uint32 acc = Accumulate();
		uint32 loEnd = 0;

		for (uint32 i = 0; i < sym_; ++i)
			loEnd += stats[i];

		rc_.EncodeFrequency(stats[sym_], loEnd, acc);

		stats[sym_] += StepSize;
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 acc = Accumulate();
		uint32 cul = rc_.GetCumulativeFreq(acc);

		uint32 idx, hiEnd;
		for (idx = 0, hiEnd = 0; (hiEnd += stats[idx]) <= cul; ++idx)
			;
		hiEnd -= stats[idx];

		rc_.UpdateFrequency(stats[idx], hiEnd, acc);
		stats[idx] += StepSize;
		return idx;
	}

private:
	static const StatType StepSize = 8;
	static const uint32 MaxAccumulatedValue = (1<<16) - 16;

	uint32 symbolCount;

	void Rescale()
	{
		for (uint32 i = 0; i < symbolCount; ++i)
			stats[i] -= stats[i] >> 1;		// no '>>=' to avoid reducing stats to 0
	}

	uint32 Accumulate()
	{
		uint32 acc = 0;
		for (uint32 i = 0; i < symbolCount; ++i)
			acc += stats[i];

		if (acc >= MaxAccumulatedValue)
		{
			Rescale();
			for (uint32 i = 0; i < symbolCount; ++i)
				acc += stats[i];
		}

		return acc;
	}

	StatType stats[MaxSymbolCount];			// can be assumed to be uint16
};

} // namespace comp

} // namespace dsrc

#endif // H_SYMBOLCODERRC
