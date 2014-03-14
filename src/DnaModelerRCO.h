/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DNAMODELRCO
#define H_DNAMODELRCO

#include "../include/dsrc/Globals.h"

#include "DnaModeler.h"
#include "Fastq.h"
#include "RangeCoder.h"
#include "BitMemory.h"
#include "SymbolCoderRC.h"

namespace dsrc
{

namespace comp
{

template <uint32 _TOrder, uint32 _TAlphabetSize>
class TDnaRCOrderModeler : public IDnaModeler
{
public:
	static const uint32 AlphabetSize = _TAlphabetSize;
	static const uint32 AlphabetBits = core::TLog2<AlphabetSize>::Value;
	static const uint32 Order = _TOrder;

	TDnaRCOrderModeler()
		:	hash(0)
	{}

	void ProcessStats(const DnaStats &stats_)
	{
		ASSERT(stats_.symbolCount <= AlphabetSize);
	}

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_)
	{
		Clear();

		RangeEncoder encoder(writer_);

		encoder.Start();
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			const fq::FastqRecord& r = records_[i];
			for (uint32 j = 0; j < r.sequenceLen; ++j)
			{
				ASSERT(r.sequence[j] < AlphabetSize);
				EncodeSymbol(encoder, r.sequence[j]);
			}
		}
		encoder.End();
	}

	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_)
	{
		Clear();

		RangeDecoder decoder(reader_);

		decoder.Start();
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			const fq::FastqRecord& r = records_[i];
			for (uint32 j = 0; j < r.sequenceLen; ++j)
			{
				r.sequence[j] = DecodeSymbol(decoder);
				ASSERT(r.sequence[j] < AlphabetSize);
			}
		}
		decoder.End();
	}


private:
	typedef uint64 HashType;
	typedef TSymbolCoderRC<AlphabetSize> Coder;
	typedef typename Coder::StatType CoderStatType;

	static const HashType HashMask = (1 << (Order * AlphabetBits)) - 1;
	static const uint32 ModelCount = 1 << (core::TLog2<AlphabetSize>::Value * Order);

	Coder coders[ModelCount];
	HashType hash;

	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		coders[GetHash()].EncodeSymbol(rc_, sym_);

		UpdateHash(sym_);
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 sym = coders[GetHash()].DecodeSymbol(rc_);

		UpdateHash(sym);

		return sym;
	}


	void Clear()
	{
		// clear hash
		hash = 0;

		// clear stats -- fill context data with initial '1' value
		CoderStatType* cd = (CoderStatType*)coders;
		std::fill(cd, cd + ModelCount * AlphabetSize, 1);
	}

	HashType GetHash()
	{
		return hash;
	}

	void UpdateHash(uint32 sym_)
	{
		hash <<= AlphabetBits;
		hash |= sym_;
		hash &= HashMask;
	}
};

} // namespace comp

} // namespace dsrc

#endif // H_DNAMODELRCO
