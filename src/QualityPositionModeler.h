/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_QUALITYPOSITIONMODELER
#define H_QUALITYPOSITIONMODELER

#include "../include/dsrc/Globals.h"

#include "QualityModeler.h"
#include "SymbolCoderRC.h"
#include "huffman.h"
#include "RecordsProcessor.h"

namespace dsrc
{

namespace comp
{

class IQualityPositionModeler : public IQualityModeler
{
public:
	IQualityPositionModeler(bool quantizedValues_)
		:	quantizedValues(quantizedValues_)
		,	symbolCount(0)
		,	minLength((uint32)-1)
		,	maxLength(0)
	{
		std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
	}

	void ProcessStats(const QualityStats &stats_)
	{
		// prepare symbol indices
		//
		symbolCount = stats_.symbolCount;
		std::copy(stats_.symbols, stats_.symbols + MaxSymbolCount, symbols);

		minLength = stats_.minLength;
		maxLength = stats_.maxLength;
	}

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_);
	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_);

protected:
	typedef HuffmanEncoder CoderType;

	static const uint32 MaxSymbolCount = QualityStats::MaxSymbolCount;
	static const byte EmptySymbol = QualityStats::EmptySymbol;

	const bool quantizedValues;

	uint32 symbolCount;
	uchar symbols[MaxSymbolCount];

	uint32 minLength;
	uint32 maxLength;

	std::vector<CoderType> positionContexts;

	virtual void EncodeRecords(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_) = 0;
	virtual void DecodeRecords(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_) = 0;
	virtual void CalculatePositionStats(const fq::FastqRecord *records_, uint32 recordsCount_, std::vector<uint32*>& stats_) = 0;

	void ComputeHuffmanContext(const fq::FastqRecord* records_, uint32 recordsCount_);

	void StoreContext(core::BitMemoryWriter& writer_);
	void ReadContext(core::BitMemoryReader& reader_);

	void StoreStatsData(core::BitMemoryWriter& writer_);
	void ReadStatsData(core::BitMemoryReader& reader_);

	void StoreSymbols(core::BitMemoryWriter& writer_);
	void ReadSymbols(core::BitMemoryReader& reader_);
};

class QualityPositionModelerPlain : public IQualityPositionModeler
{
public:
	QualityPositionModelerPlain(bool quantizedValues_)
		:	IQualityPositionModeler(quantizedValues_)
	{}

private:
	void CalculatePositionStats(const fq::FastqRecord *records_, uint32 recordsCount_, std::vector<uint32 *> &stats_);
	void EncodeRecords(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_);
	void DecodeRecords(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_);
};


class QualityPositionModelerTruncated : public IQualityPositionModeler
{
public:
	QualityPositionModelerTruncated(bool quantizedValues_)
		:	IQualityPositionModeler(quantizedValues_)
	{}

private:
	void CalculatePositionStats(const fq::FastqRecord *records_, uint32 recordsCount_, std::vector<uint32 *> &stats_);
	void EncodeRecords(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_);
	void DecodeRecords(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_);
};

} // namespace comp

} // namespace dsrc

#endif // H_QUALITYPOSITIONMODELER
