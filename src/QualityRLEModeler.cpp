/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "QualityRLEModeler.h"
#include "Fastq.h"

namespace dsrc
{

namespace comp
{

using namespace core;
using namespace fq;

// QualityRLEModeler
//
void QualityRLEModeler::ProcessStats(const QualityStats &stats_)
{
	rawRunLength = stats_.rawLength;
}

void QualityRLEModeler::StoreStatsData(BitMemoryWriter &writer_)
{
	writer_.PutWord(runLength);
}

void QualityRLEModeler::StoreRunsSymbols(BitMemoryWriter &writer_)
{
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		writer_.PutBit(qSymbols[i] != EmptySymbol);
	}
	for (uint32 i = 0; i < MaxLengthSymbolsCount; ++i)
	{
		writer_.PutBit(lSymbols[i] != EmptySymbol);
	}
}

void QualityRLEModeler::ReadStatsData(BitMemoryReader &reader_)
{
	runLength = reader_.GetWord();
}

void QualityRLEModeler::ReadRunsSymbols(BitMemoryReader &reader_)
{
	// read symbols
	//
	qSymbolCount = 0;
	lSymbolCount = 0;
	std::fill(qSymbols, qSymbols + MaxSymbolCount, +EmptySymbol);
	std::fill(lSymbols, lSymbols + MaxLengthSymbolsCount, +EmptySymbol);
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		if (reader_.GetBit())
		{
			qSymbols[qSymbolCount++] = i;
		}
	}
	for (uint32 i = 0; i < MaxLengthSymbolsCount; ++i)
	{
		if (reader_.GetBit())
		{
			lSymbols[lSymbolCount++] = i;
		}
	}
	reader_.FlushInputWordBuffer();
	ASSERT(qSymbolCount > 0);
	ASSERT(lSymbolCount > 0);
}

void QualityRLEModeler::DecodeRecords(FastqRecord *records_, uint32 recordsCount_)
{
	uchar* symRun = symbolRun.Pointer();
	uchar* lenRun = lengthsRun.Pointer();

	uint32 curLen = 0;
	uchar curQua = 0;
	uint32 curIdx = 0;
	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		FastqRecord& r = records_[i];
		uint32 nCount = 0;

		for (uint32 j = 0; j < r.qualityLen; ++j)
		{
			if (curLen == 0)
			{
				ASSERT(curIdx < runLength);

				ASSERT(symRun[curIdx] <= MaxQualitySymbol);
				curQua = symRun[curIdx];

				ASSERT(lenRun[curIdx] <= MaxLengthSymbol);
				curLen = lenRun[curIdx] + 1;
				curIdx++;
			}
			r.quality[j] = curQua;
			--curLen;

			if (quantizedValues)
				nCount += (uint32)(curQua == 0);
			else
				nCount += (uint32)(curQua >= 128);
		}

		// remember to update sequence length !
		r.sequenceLen = r.qualityLen - nCount;
	}
}


// QualityRLEModeler
//
void QualityRLEModeler::Encode(BitMemoryWriter &writer_, const FastqRecord *records_, uint32 recordsCount_)
{
	EncodeRecords(records_, recordsCount_);

	CalculateSymbolIndices();

	ComputeHuffmanContext();

	writer_.FlushPartialWordBuffer();

	StoreStatsData(writer_);

	StoreRunsSymbols(writer_);

	StoreHuffmanContext(writer_);

	EncodeRuns(writer_);

	writer_.FlushPartialWordBuffer();
}

void QualityRLEModeler::EncodeRecords(const FastqRecord *records_, uint32 recordsCount_)
{
	// prepare buffers
	//
	if (symbolRun.Size() < rawRunLength)
		symbolRun.Extend(rawRunLength);
	if (lengthsRun.Size() < rawRunLength)
		lengthsRun.Extend(rawRunLength);

	uchar* symRun = symbolRun.Pointer();
	uchar* lenRun = lengthsRun.Pointer();

	uint32* qF = (uint32*)qFreqs.Pointer();
	uint32* lF = (uint32*)lFreqs.Pointer();


	// clear data
	//
	runLength = 0;

	ASSERT(qFreqs.Size() >= MaxSymbolCount * sizeof(uint32));
	ASSERT(lFreqs.Size() >= MaxLengthSymbolsCount * sizeof(uint32));
	std::fill(qF, qF + MaxSymbolCount, 0);
	std::fill(lF, lF + MaxLengthSymbolsCount, 0);

	uchar prevSym = EmptySymbol;
	uchar curLen = 0;


	// calculate runs and check available symbols
	//
	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		const FastqRecord& r = records_[i];
		for (uint32 j = 0; j < r.qualityLen; ++j)
		{
			if (r.quality[j] == prevSym && curLen < MaxLengthSymbol)
			{
				curLen++;
			}
			else
			{
				if (prevSym != EmptySymbol)
				{
					symRun[runLength  ] = prevSym;
					lenRun[runLength++] = curLen;

					qF[prevSym]++;
					lF[curLen]++;
				}

				curLen = 0;
				prevSym = r.quality[j];
			}
		}
	}
	ASSERT(curLen < MaxLengthSymbol);

	symRun[runLength  ] = prevSym;
	lenRun[runLength++] = curLen;

	qF[prevSym]++;
	lF[curLen]++;
}

void QualityRLEModeler::CalculateSymbolIndices()
{
	std::fill(lSymbols, lSymbols + MaxLengthSymbolsCount, +EmptySymbol);
	std::fill(qSymbols, qSymbols + MaxSymbolCount, +EmptySymbol);

	uint32* qF = (uint32*)qFreqs.Pointer();
	uint32* lF = (uint32*)lFreqs.Pointer();

	qSymbolCount = 0;
	lSymbolCount = 0;
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		if (qF[i] != 0)
		{
			qSymbols[i] = qSymbolCount++;
		}

		if (lF[i] != 0)
		{
			lSymbols[i] = lSymbolCount++;
		}
	}
	ASSERT(lSymbolCount > 0);
	ASSERT(qSymbolCount > 0);
}

void QualityRLEModeler::ComputeHuffmanContext()
{
	if (qSymbolCount == 1)
		return;

	// prepare buffers
	//
	{
		uint32 qfs = qSymbolCount * qSymbolCount * sizeof(uint32);
		uint32 lfs = qSymbolCount * lSymbolCount * sizeof(uint32);

		if (qFreqs.Size() < qfs)
			qFreqs.Extend(qfs);
		if (lFreqs.Size() < lfs)
			lFreqs.Extend(lfs);

		std::fill(qFreqs.Pointer(), qFreqs.Pointer() + qfs, 0);
		std::fill(lFreqs.Pointer(), lFreqs.Pointer() + lfs, 0);
	}

	std::vector<uint32*> qF;
	qF.resize(qSymbolCount);
	for (uint32 i = 0; i < qSymbolCount; ++i)
		qF[i] = (uint32*)qFreqs.Pointer() + i * qSymbolCount;

	std::vector<uint32*> lF;
	lF.resize(qSymbolCount);
	for (uint32 i = 0; i < qSymbolCount; ++i)
		lF[i] = (uint32*)lFreqs.Pointer() + i * lSymbolCount;

	uchar* qRun = symbolRun.Pointer();
	uchar* lRun = lengthsRun.Pointer();


	// compute stats
	//
	uchar prev = 0;
	for (uint32 i = 0; i < runLength; ++i)
	{
		ASSERT(prev < qSymbolCount);

		uint32 q = qSymbols[qRun[i]];
		ASSERT(q < qSymbolCount);

		uint32 l = lSymbols[lRun[i]];
		ASSERT(l < lSymbolCount);

		qF[prev][q]++;
		lF[q][l]++;
		prev = q;
	}


	// compute huffman context
	//
	qContexts.clear();
	qContexts.resize(qSymbolCount);
	lContexts.clear();
	lContexts.resize(qSymbolCount);

	for (uint32 i = 0; i < qSymbolCount; ++i)
	{
		qContexts[i].Restart(qSymbolCount);
		for (uint32 j = 0; j < qSymbolCount; ++j)
		{
			qContexts[i].Insert(qF[i][j]);
		}
		qContexts[i].Complete();

		lContexts[i].Restart(lSymbolCount);
		for (uint32 j = 0; j < lSymbolCount; ++j)
		{
			lContexts[i].Insert(lF[i][j]);
		}
		lContexts[i].Complete();
	}

}

void QualityRLEModeler::StoreHuffmanContext(BitMemoryWriter &writer_)
{
	if (qSymbolCount == 1)
		return;

	for (uint32 i = 0; i < qSymbolCount; ++i)
	{
		qContexts[i].StoreTree(writer_);
		lContexts[i].StoreTree(writer_);
	}
}

void QualityRLEModeler::EncodeRuns(BitMemoryWriter &writer_)
{
	if (qSymbolCount > 1)
	{
		// get huffman codes
		//
		std::vector<HuffmanEncoder::Code const*> qCodes, lCodes;
		qCodes.resize(qSymbolCount);
		lCodes.resize(qSymbolCount);

		for (uint32 i = 0; i < qSymbolCount; ++i)
		{
			qCodes[i] = qContexts[i].GetCodes();
			lCodes[i] = lContexts[i].GetCodes();
		}

		const uchar* qRun = symbolRun.Pointer();
		const uchar* lRun = lengthsRun.Pointer();

		// encode runs
		//
		uchar prev = 0;
		for (uint32 i = 0; i < runLength; ++i)
		{
			uchar q = qSymbols[qRun[i]];
			uchar l = lSymbols[lRun[i]];

			const HuffmanEncoder::Code& qCode = qCodes[prev][q];
			const HuffmanEncoder::Code& lCode = lCodes[q][l];

			writer_.PutBits(qCode.code, qCode.len);
			writer_.PutBits(lCode.code, lCode.len);

			prev = q;
		}
	}
	else
	{
		ASSERT(lSymbolCount <= 2);

		if (lSymbolCount > 1)
		{
			writer_.FlushPartialWordBuffer();

			const uchar* lRun = lengthsRun.Pointer();
			uchar begin = lSymbols[lRun[0]];
			writer_.PutByte(begin);
		}
	}
}


void QualityRLEModeler::Decode(BitMemoryReader &reader_, FastqRecord *records_, uint32 recordsCount_)
{
	ReadStatsData(reader_);

	ReadRunsSymbols(reader_);

	ReadHuffmanContext(reader_);

	reader_.FlushInputWordBuffer();

	DecodeRuns(reader_);

	DecodeRecords(records_, recordsCount_);

	reader_.FlushInputWordBuffer();
}

void QualityRLEModeler::ReadHuffmanContext(BitMemoryReader &reader_)
{
	if (qSymbolCount == 1)
		return;

	qContexts.clear();
	qContexts.resize(qSymbolCount);
	lContexts.clear();
	lContexts.resize(qSymbolCount);

	for (uint32 i = 0; i < qSymbolCount; ++i)
	{
		qContexts[i].LoadTree(reader_);
		lContexts[i].LoadTree(reader_);
	}
}

void QualityRLEModeler::DecodeRuns(BitMemoryReader &reader_)
{
	// prepare buffers
	//
	if (symbolRun.Size() < runLength)
		symbolRun.Extend(runLength);
	if (lengthsRun.Size() < runLength)
		lengthsRun.Extend(runLength);

	uchar* qRun = symbolRun.Pointer();
	uchar* lRun = lengthsRun.Pointer();

	if (qSymbolCount > 1)
	{
		uint32 prev = 0;
		for (uint32 i = 0; i < runLength; ++i)
		{
			// decode quality
			uint32 bit = reader_.GetBits(qContexts[prev].GetMinLen());
			int32 idx = qContexts[prev].DecodeFast(bit);

			while (idx < 0)
			{
				bit = reader_.GetBit();
				idx = qContexts[prev].Decode(bit);
			};

			ASSERT((uint32)idx < qSymbolCount);
			ASSERT(qSymbols[idx] != EmptySymbol);
			qRun[i] = qSymbols[idx];
			prev = idx;

			// decode length
			bit = reader_.GetBits(lContexts[prev].GetMinLen());
			idx = lContexts[prev].DecodeFast(bit);

			while (idx < 0)
			{
				bit = reader_.GetBit();
				idx = lContexts[prev].Decode(bit);
			};

			ASSERT((uint32)idx < lSymbolCount);
			ASSERT(lSymbols[idx] != EmptySymbol);
			lRun[i] = lSymbols[idx];
		}
	}
	else
	{
		uchar qSym = qSymbols[0];
		uchar lBegin = 0;
		uchar lEnd = 0;
		if (lSymbolCount > 1)
		{
			reader_.FlushInputWordBuffer();

			lBegin = reader_.GetByte();
			lBegin = lSymbols[lBegin];
			lEnd = lSymbols[0];
			if (lEnd == lBegin)
				lEnd = lSymbols[1];
			ASSERT(lEnd != lBegin);
		}
		else
		{
			lBegin = lSymbols[0];
			lEnd = lBegin;
		}

		ASSERT(lBegin != EmptySymbol);
		ASSERT(lEnd != EmptySymbol);
		std::fill(qRun, qRun + runLength, qSym);
		std::fill(lRun, lRun + runLength, lBegin);
		lRun[runLength - 1] = lEnd;
	}
}

} // namespace comp

} // namespace dsrc
