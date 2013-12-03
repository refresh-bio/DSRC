/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include "QualityPositionModeler.h"
#include "Fastq.h"

namespace dsrc
{

namespace comp
{

using namespace fq;
using namespace core;

// IQualityPositionModeler
//
void IQualityPositionModeler::StoreStatsData(BitMemoryWriter &writer_)
{
	writer_.PutWord(maxLength);
}

void IQualityPositionModeler::StoreSymbols(BitMemoryWriter &writer_)
{
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		writer_.PutBit(symbols[i] != (uchar)EmptySymbol);
	}
}

void IQualityPositionModeler::ReadStatsData(BitMemoryReader &reader_)
{
	maxLength = reader_.GetWord();
	ASSERT(maxLength > 0);
}

void IQualityPositionModeler::ReadSymbols(BitMemoryReader &reader_)
{
	symbolCount = 0;
	std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		if (reader_.GetBit())
		{
			symbols[symbolCount++] = i;
		}
	}
	ASSERT(symbolCount > 0);
}

void IQualityPositionModeler::Encode(BitMemoryWriter& writer_, const FastqRecord* records_, uint32 recordsCount_)
{
	ComputeHuffmanContext(records_, recordsCount_);

	writer_.FlushPartialWordBuffer();

	StoreStatsData(writer_);

	StoreSymbols(writer_);

	StoreContext(writer_);

	EncodeRecords(writer_, records_, recordsCount_);

	writer_.FlushPartialWordBuffer();
}

void IQualityPositionModeler::StoreContext(BitMemoryWriter &writer_)
{
	for (uint32 i = 0; i < maxLength; ++i)
	{
		positionContexts[i].StoreTree(writer_);
	}
}

void IQualityPositionModeler::Decode(BitMemoryReader &reader_, FastqRecord *records_, uint32 recordsCount_)
{
	reader_.FlushInputWordBuffer();

	ReadStatsData(reader_);

	ReadSymbols(reader_);

	ReadContext(reader_);

	DecodeRecords(reader_, records_, recordsCount_);

	reader_.FlushInputWordBuffer();
}

void IQualityPositionModeler::ReadContext(BitMemoryReader &reader_)
{
	positionContexts.clear();
	positionContexts.resize(maxLength);
	for (uint32 i = 0; i < maxLength; ++i)
	{
		positionContexts[i].LoadTree(reader_);
	}
}

void IQualityPositionModeler::ComputeHuffmanContext(const FastqRecord *records_, uint32 recordsCount_)
{
	// initialize position stats
	//
	Buffer statsBuffer(maxLength * symbolCount * sizeof(uint32));
	std::fill(statsBuffer.Pointer(), statsBuffer.Pointer() + statsBuffer.Size(), 0);

	std::vector<uint32*> positionStats;
	positionStats.resize(maxLength);
	for (uint32 i = 0; i < maxLength; ++i)
	{
		positionStats[i] = (uint32*)(statsBuffer.Pointer() + i * symbolCount * sizeof(uint32));
	}

	// calculate position stats
	//
	CalculatePositionStats(records_, recordsCount_, positionStats);

	// initialize context
	//
	positionContexts.clear();				// there's a nasty bug with HuffmanEncoder...
	positionContexts.resize(maxLength);
	for (uint32 i = 0; i < maxLength; ++i)
	{
		positionContexts[i].Restart(symbolCount);
		for (uint32 j = 0; j < symbolCount; ++j)
		{
			positionContexts[i].Insert(positionStats[i][j]);
		}
		positionContexts[i].Complete();
	}
}


// QualityPositionModelerPlain
//
void QualityPositionModelerPlain::CalculatePositionStats(const FastqRecord *records_, uint32 recordsCount_, std::vector<uint32 *> &stats_)
{
	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		const FastqRecord& r = records_[i];

		ASSERT(stats_.size() >= r.qualityLen);

		for (uint32 j = 0; j < r.qualityLen; ++j)
		{
			ASSERT(symbols[r.quality[j]] != EmptySymbol);

			stats_[j][symbols[r.quality[j]]]++;
		}
	}
}

void QualityPositionModelerPlain::EncodeRecords(BitMemoryWriter &writer_, const FastqRecord *records_, uint32 recordsCount_)
{
	// get coding symbols and store trees
	//
	ASSERT(positionContexts.size() == maxLength);
	std::vector<HuffmanEncoder::Code*> hufCodes;
	hufCodes.resize(maxLength);

	for (uint32 i = 0; i < maxLength; ++i)
	{
		hufCodes[i] = (HuffmanEncoder::Code*)positionContexts[i].GetCodes();
	}

	// store records
	//
	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		const FastqRecord& r = records_[i];
		for (uint32 j = 0; j < r.qualityLen; ++j)
		{
			ASSERT(symbols[r.quality[j]] != EmptySymbol);

			int32 qua = symbols[r.quality[j]];
			ASSERT(qua < (int32)symbolCount);
			writer_.PutBits(hufCodes[j][qua].code, hufCodes[j][qua].len);
		}
	}
}

void QualityPositionModelerPlain::DecodeRecords(BitMemoryReader &reader_, FastqRecord *records_, uint32 recordsCount_)
{

	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		FastqRecord& r = records_[i];

		uint32 nCount = 0;
		for (uint32 j = 0; j < r.qualityLen; ++j)
		{
			uint32 bit = reader_.GetBits(positionContexts[j].GetMinLen());
			int32 idx = positionContexts[j].DecodeFast(bit);

			while (idx < 0)
			{
				bit = reader_.GetBit();
				idx = positionContexts[j].Decode(bit);
			};

			r.quality[j] = symbols[idx];

			if (quantizedValues)
				nCount += (uint32)(r.quality[j] == 0);
			else
				nCount += (uint32)(r.quality[j] >= 128);
		}

		r.sequenceLen = r.qualityLen - nCount;
	}
}


// IQualityPositionTruncated
//
void QualityPositionModelerTruncated::CalculatePositionStats(const FastqRecord *records_, uint32 recordsCount_, std::vector<uint32 *> &stats_)
{
	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		const FastqRecord& r = records_[i];

		ASSERT(stats_.size() >= r.truncatedLen);

		for (uint32 j = 0; j < r.truncatedLen; ++j)
		{
			ASSERT(symbols[r.quality[j]] != EmptySymbol);

			stats_[j][symbols[r.quality[j]]]++;
		}
	}
}

void QualityPositionModelerTruncated::EncodeRecords(BitMemoryWriter &writer_, const FastqRecord *records_, uint32 recordsCount_)
{
	// get coding symbols and store trees
	//
	ASSERT(positionContexts.size() == maxLength);
	std::vector<HuffmanEncoder::Code*> hufCodes;
	hufCodes.resize(maxLength);

	for (uint32 i = 0; i < maxLength; ++i)
	{
		hufCodes[i] = (HuffmanEncoder::Code*)positionContexts[i].GetCodes();
	}


	// store records
	//
	const bool variableLength = minLength != maxLength;
	const uint32 maxBitLength = core::bit_length(maxLength);

	writer_.PutBit(variableLength);
	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		const FastqRecord& r = records_[i];

		// store info about truncation
		//
		ASSERT(r.truncatedLen <= r.qualityLen);
		writer_.PutBit(r.qualityLen != r.truncatedLen);
		if (r.qualityLen != r.truncatedLen)
		{
			//uint32 diff = r.qualityLen - r.truncatedLen;
			uint32 bitLen = (variableLength) ? core::bit_length(r.qualityLen) : maxBitLength;
			writer_.PutBits(r.truncatedLen, bitLen);
		}

		// store quality
		//
		ASSERT(r.qualityLen >= r.truncatedLen);
		for (uint32 j = 0; j < r.truncatedLen; ++j)
		{
			ASSERT(symbols[r.quality[j]] != EmptySymbol);

			int32 qua = symbols[r.quality[j]];
			ASSERT(qua < (int32)symbolCount);
			writer_.PutBits(hufCodes[j][qua].code, hufCodes[j][qua].len);
		}
	}
}

void QualityPositionModelerTruncated::DecodeRecords(BitMemoryReader &reader_, FastqRecord *records_, uint32 recordsCount_)
{
	const uint32 maxBitLength = core::bit_length(maxLength);
	const bool variableLength = reader_.GetBit();

	const byte HashSymbol = quantizedValues ? +HashSymbolQuantized : +HashSymbolNormal;

	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		FastqRecord& r = records_[i];

		uint32 thLen = r.qualityLen;
		if (reader_.GetBit())
		{
			uint32 bitLen = (variableLength) ? core::bit_length(r.qualityLen) : maxBitLength;
			thLen = reader_.GetBits(bitLen);
		}
		ASSERT(thLen <= r.qualityLen);

		uint32 nCount = 0;
		for (uint32 j = 0; j < thLen; ++j)
		{
			uint32 bit = reader_.GetBits(positionContexts[j].GetMinLen());
			int32 idx = positionContexts[j].DecodeFast(bit);

			while (idx < 0)
			{
				bit = reader_.GetBit();
				idx = positionContexts[j].Decode(bit);
			};

			ASSERT(idx < (int32)MaxSymbolCount);
			r.quality[j] = symbols[idx];

			if (quantizedValues)
				nCount += (uint32)(r.quality[j] == 0);
			else
				nCount += (uint32)(r.quality[j] >= 128);
		}

		for (uint32 j = thLen; j < r.qualityLen; ++j)
		{
			r.quality[j] = HashSymbol;			// hash symbol
		}

		r.sequenceLen = r.qualityLen - nCount;
	}
}

} // namespace comp

} // namespace dsrc
