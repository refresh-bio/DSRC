/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_QUALITYENCODER
#define H_QUALITYENCODER

#include "../include/dsrc/Globals.h"

#include "RangeCoder.h"
#include "SymbolCoderRC.h"

namespace dsrc
{

namespace comp
{

template <uint32 _TSymbolCount, uint32 _TSymbolOrder, uint32 _TTotalOrder = _TSymbolOrder>
class TQualityModelBase
{
public:
	static const uint32 SymbolCount = _TSymbolCount;
	static const uint32 SymbolOrder = _TSymbolOrder;
	static const uint32 TotalOrder = _TTotalOrder;

	TQualityModelBase()
		:	models(NULL)
		,	hash(0)
		,	symBuffer(0)
	{
		models = new Coder[ModelCount];
	}

	~TQualityModelBase()
	{
		delete models;
	}

	void Clear()
	{
		hash = 0;
		symBuffer = 0;
		CoderStatType* cd = (CoderStatType*)models;
		std::fill(cd, cd + ModelCount * SymbolCount, 1);
	}

protected:
	typedef uint64 THash;

	static const uint32 AlphabetSize = SymbolCount;
	static const uint32	AlphabetBits = core::TLog2<AlphabetSize>::Value;
	static const uint32	ModelCount = 1 << (AlphabetBits * TotalOrder);

	static const THash HashMask = (1ULL << (TotalOrder * AlphabetBits)) - 1ULL;
	static const THash SymbolMask = (1ULL << AlphabetBits) - 1ULL;

	static const uint64 BitsLo = (SymbolOrder/2) * AlphabetBits;
	static const uint64 BitsHi = (SymbolOrder/2 + 1) * AlphabetBits;
	static const THash SymbolSwapMask = core::TBitMask<BitsLo>::Value | ~core::TBitMask<BitsHi>::Value;

	static const THash SymbolHashMask = (1ULL << (SymbolOrder * AlphabetBits)) - 1ULL;
	static const THash SymbolContextBits = AlphabetBits * SymbolOrder;

	typedef TSymbolCoderRC<AlphabetSize> Coder;
	typedef typename Coder::StatType CoderStatType;

	Coder* models;
	THash hash;
	THash symBuffer;

	void UpdateHash(uint32 sym_)
	{
		hash <<= AlphabetBits;

		uint64 nextBuf = ((hash >> BitsLo) & SymbolMask);
		uint64 swp = (nextBuf + symBuffer) / 2;

		hash &= SymbolSwapMask;
		hash |= (swp << BitsLo);
		hash |= sym_;

		symBuffer = nextBuf;
	}

	THash GetHash() const
	{
		return hash & SymbolHashMask;
	}
};


template <uint32 _TSymbolCount, uint32 _TOrder>
class TQualityModel : public TQualityModelBase<_TSymbolCount, _TOrder, _TOrder>
{
public:
	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		Super::models[Super::GetHash()].EncodeSymbol(rc_, sym_);

		Super::UpdateHash(sym_);
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 sym = Super::models[Super::GetHash()].DecodeSymbol(rc_);

		Super::UpdateHash(sym);
		return sym;
	}

private:
	typedef TQualityModelBase<_TSymbolCount, _TOrder, _TOrder> Super;
};


template <uint32 _TSymbolCount, uint32 _TOrder>
class TQualityModelExt : public TQualityModelBase<_TSymbolCount, _TOrder, _TOrder + 1>
{
public:
	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_, uint32 ctx0_)
	{
		ASSERT(ctx0_ < Super::AlphabetSize);
		ASSERT(sym_ < Super::AlphabetSize);

		uint32 h = (Super::GetHash() << Super::AlphabetBits) | ctx0_;

		Super::models[h].EncodeSymbol(rc_, sym_);

		Super::UpdateHash(sym_);
	}

	uint32 DecodeSymbol(RangeDecoder& rc_, uint32 ctx0_)
	{
		ASSERT(ctx0_ < Super::AlphabetSize);

		uint32 h = (Super::GetHash()  << Super::AlphabetBits) | ctx0_;
		uint32 sym = Super::models[h].DecodeSymbol(rc_);

		Super::UpdateHash(sym);
		return sym;
	}

private:
	typedef TQualityModelBase<_TSymbolCount, _TOrder, _TOrder + 1> Super;
};


struct SpecialSymbolNormalHandlerStrategy
{
	static bool IsValid(byte symbol_)
	{
		return symbol_ >= 128;
	}
};

struct SpecialSymbolLossyHandlerStrategy
{
	static bool IsValid(byte symbol_)
	{
		return symbol_ == 0;
	}
};


template <class _TModel, class _TSpecialSymbolHandler>
class TNormalQualityEncoder
{
public:
	typedef _TModel Model;
	typedef _TSpecialSymbolHandler SymbolHandler;

	void ProcessStats(const QualityStats& stats_)
	{
		ASSERT(stats_.symbolCount < SymbolCount);
	}

	void Encode(const fq::FastqRecord& rec_, Model& model_, RangeEncoder& coder_)
	{
		for (uint32 j = 0; j < rec_.qualityLen; ++j)
		{
			uint32 ctx0 = rec_.quality[j];
			ASSERT(ctx0 < SymbolCount);
			ASSERT(ctx0 != EmptySymbol);

			model_.EncodeSymbol(coder_, ctx0);
		}
	}

	void Decode(fq::FastqRecord& rec_, Model& model_, RangeDecoder& coder_)
	{
		uint32 nCount = 0;

		for (uint32 j = 0; j < rec_.qualityLen; ++j)
		{
			uint32 c = model_.DecodeSymbol(coder_);
			ASSERT(c < SymbolCount);

			rec_.quality[j] = c;
			nCount += SymbolHandler::IsValid(c);
			//nCount += (uint32)(c == 0);
		}

		rec_.sequenceLen = rec_.qualityLen - nCount;
	}

	void Store(core::BitMemoryWriter& )
	{}

	void Read(core::BitMemoryReader& )
	{}

private:
	static const uint32 SymbolCount = Model::SymbolCount;
	static const byte EmptySymbol = 255;
};


template <class _TModel, class _TSpecialSymbolHandler>
class TPositionalQualityEncoder
{
public:
	typedef _TModel Model;
	typedef _TSpecialSymbolHandler SymbolHandler;

	void ProcessStats(const QualityStats& stats_)
	{
		ASSERT(stats_.symbolCount < SymbolCount);
	}

	void Encode(const fq::FastqRecord& rec_, Model& model_, RangeEncoder& coder_)
	{
		for (uint32 j = 0; j < rec_.qualityLen; ++j)
		{
			uint32 ctx0 = rec_.quality[j];
			ASSERT(ctx0 < SymbolCount);
			ASSERT(ctx0 != EmptySymbol);

			uint32 pctx = (j * SymbolCount / rec_.qualityLen);
			ASSERT(pctx < SymbolCount);
			model_.EncodeSymbol(coder_, ctx0, pctx);
		}
	}

	void Decode(fq::FastqRecord& rec_, Model& model_, RangeDecoder& coder_)
	{
		uint32 nCount = 0;

		for (uint32 j = 0; j < rec_.qualityLen; ++j)
		{
			uint32 pctx = (j * SymbolCount / rec_.qualityLen);
			ASSERT(pctx < SymbolCount);

			uint32 c = model_.DecodeSymbol(coder_, pctx);
			ASSERT(c < SymbolCount);

			rec_.quality[j] = c;
			nCount += SymbolHandler::IsValid(c);
			//nCount += (uint32)(c == 0);
		}

		rec_.sequenceLen = rec_.qualityLen - nCount;
	}

	void Store(core::BitMemoryWriter& )
	{}

	void Read(core::BitMemoryReader& )
	{}

private:
	static const uint32 SymbolCount = Model::SymbolCount;
	static const byte EmptySymbol = 255;
};

template <class _TModel, class _TSpecialSymbolHandler, uint32 _TSymbolRescale>
class TTranslationalQualityEncoder
{
public:
	typedef _TModel Model;
	typedef _TSpecialSymbolHandler SymbolHandler;

	TTranslationalQualityEncoder()
	{
		std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
	}

	void ProcessStats(const QualityStats& stats_)
	{
		ASSERT(stats_.symbolCount <= SymbolCount);
		std::copy(stats_.symbols, stats_.symbols + MaxSymbolCount, symbols);
	}

	void Encode(const fq::FastqRecord& rec_, Model& model_, RangeEncoder& coder_)
	{
		for (uint32 j = 0; j < rec_.qualityLen; ++j)
		{
			uint32 ctx0 = symbols[rec_.quality[j]];
			ASSERT(ctx0 < SymbolCount);
			ASSERT(ctx0 != EmptySymbol);

			uint32 pctx = (j * SymbolRescale / rec_.qualityLen);
			ASSERT(pctx < SymbolCount);
			model_.EncodeSymbol(coder_, ctx0, pctx);
		}
	}

	void Decode(fq::FastqRecord& rec_, Model& model_, RangeDecoder& coder_)
	{
		uint32 nCount = 0;

		for (uint32 j = 0; j < rec_.qualityLen; ++j)
		{
			uint32 pctx = (j * SymbolRescale/ rec_.qualityLen);
			ASSERT(pctx < SymbolCount);

			uint32 c = model_.DecodeSymbol(coder_, pctx);
			ASSERT(c < SymbolCount);

			//ASSERT(c < symbolCount);
			rec_.quality[j] = symbols[c];
			nCount += SymbolHandler::IsValid(symbols[c]);		}

		rec_.sequenceLen = rec_.qualityLen - nCount;
	}

	void Store(core::BitMemoryWriter& writer_)
	{
		writer_.FlushPartialWordBuffer();

		// this method is more sufficient where SymbolCount >= 64, as we store 256 bits (32B)
		for (uint32 i = 0; i < MaxSymbolCount; ++i)
		{
			writer_.PutBit(symbols[i] != EmptySymbol);
		}
		writer_.FlushFullWordBuffer();
	}

	void Read(core::BitMemoryReader& reader_)
	{
		std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);

		reader_.FlushInputWordBuffer();

		// this method is more sufficient where SymbolCount >= 64, as we store 256 bits (32B)
		uint32 symCount = 0;
		for (uint32 i = 0; i < MaxSymbolCount; ++i)
		{
			if (reader_.GetBit() != 0)
				symbols[symCount++] = i;
		}
		reader_.FlushInputWordBuffer();
	}

private:
	static const uint32 MaxSymbolCount = 256;
	static const uint32 SymbolCount = Model::SymbolCount;
	static const uint32 SymbolRescale = _TSymbolRescale;
	static const byte EmptySymbol = 255;

	byte symbols[MaxSymbolCount];
};

} // namespace comp

} // namespace dsrc

#endif // H_QUALITYENCODER
