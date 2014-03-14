/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_QUALITYRLEMODELER
#define H_QUALITYRLEMODELER

#include "../include/dsrc/Globals.h"

#include "Fastq.h"
#include "QualityModeler.h"
#include "SymbolCoderRC.h"

#include "huffman.h"

namespace dsrc
{

namespace comp
{

class QualityRLEModeler : public IQualityModeler
{
public:
	QualityRLEModeler(bool quantizedValues_)
		:	quantizedValues(quantizedValues_)
		,	qSymbolCount(0)
		,	lSymbolCount(0)
		,	rawRunLength(0)
		,	runLength(0)
		,	symbolRun(8)
		,	lengthsRun(8)
		,	qFreqs(MaxSymbolCount * sizeof(uint32))
		,	lFreqs(MaxLengthSymbolsCount * sizeof(uint32))
	{
		ClearSymbols();
	}

	void ProcessStats(const QualityStats& stats_);

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_);
	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_);

protected:
	typedef HuffmanEncoder CoderType;

	static const uint32 MaxSymbolCount = QualityStats::MaxSymbolCount;
	static const byte EmptySymbol = QualityStats::EmptySymbol;

	static const uint32 MaxLengthSymbolsCount = 256;
	static const uint32 MaxQualitySymbol = 254;
	static const uint32 MaxLengthSymbol = 254;

	const bool quantizedValues;

	uint32 qSymbolCount;
	uchar qSymbols[MaxSymbolCount];
	uint32 lSymbolCount;
	uchar lSymbols[MaxLengthSymbolsCount];

	uint32 rawRunLength;
	uint32 runLength;
	core::Buffer symbolRun;
	core::Buffer lengthsRun;

	std::vector<HuffmanEncoder> lContexts;
	std::vector<HuffmanEncoder> qContexts;

	core::Buffer qFreqs;
	core::Buffer lFreqs;

	void StoreStatsData(core::BitMemoryWriter& writer_);
	void ReadStatsData(core::BitMemoryReader& reader_);

	void StoreRunsSymbols(core::BitMemoryWriter& writer_);
	void ReadRunsSymbols(core::BitMemoryReader& reader_);

	void EncodeRecords(const fq::FastqRecord* records_, uint32 recordsCount_);
	void DecodeRecords(fq::FastqRecord* records_, uint32 recordsCount_);

	void ComputeHuffmanContext();
	void CalculateSymbolIndices();

	void StoreHuffmanContext(core::BitMemoryWriter& writer_);
	void ReadHuffmanContext(core::BitMemoryReader& reader_);

	void EncodeRuns(core::BitMemoryWriter& writer_);
	void DecodeRuns(core::BitMemoryReader& reader_);

	void ClearSymbols()
	{
		std::fill(qSymbols, qSymbols + MaxSymbolCount, +EmptySymbol);
		std::fill(lSymbols, lSymbols + MaxLengthSymbolsCount, +EmptySymbol);
	}
};

} // namespace comp

} // namespace dsrc

#endif // H_QUALITYRLEMODELER
