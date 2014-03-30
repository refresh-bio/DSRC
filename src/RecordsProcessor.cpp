/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "RecordsProcessor.h"
#include "Common.h"
#include "BitMemory.h"
#include "Crc32.h"

namespace dsrc
{

namespace comp
{

using namespace fq;

// IRecordsProcessor
//
void IRecordsProcessor::ProcessRecordFromColorSpace(FastqRecord& rec_)
{
	// TODO: move to class + create LUT
	const char deltas[] = {'N', 'N', 'A', 'C', 'G', 'T',
						   'N', 'N', 'C', 'A', 'T', 'G',
						   'N', 'N', 'G', 'T', 'A', 'C',
						   'N', 'N', 'T', 'G', 'C', 'A'};

	const char* deltaA = deltas;
	const char* deltaC = deltas + 6;
	const char* deltaG = deltas + 12;
	const char* deltaT = deltas + 18;
	const char* lastMatrix = deltaA;

	uchar symbol = rec_.sequence[0];

	for (uint32 k = 1; k < rec_.sequenceLen; ++k)
	{
		// TODO: use LUT to eliminate switch
		switch (symbol)
		{
			case 'A': lastMatrix = deltaA; break;
			case 'C': lastMatrix = deltaC; break;
			case 'G': lastMatrix = deltaG; break;
			case 'T': lastMatrix = deltaT; break;
			case 'N': // symbol undefined/error, use previous valid symbol matrix
			default : break;
		}
		symbol = lastMatrix[rec_.sequence[k] - '.'];
		ASSERT(symbol =='A' || symbol == 'C' || symbol == 'G' || symbol == 'T' || symbol == 'N');

		rec_.sequence[k] = symbol;
	}
}

void IRecordsProcessor::ProcessRecordToColorSpace(FastqRecord& rec_, bool useConstDelta_, uchar seqStart_, uchar quaStart_)
{
	const char deltas[] = {'N', 'N', 'A', 'C', 'G', 'T',
						  'N', 'N', 'C', 'A', 'T', 'G',
						  'N', 'N', 'G', 'T', 'A', 'C',
						  'N', 'N', 'T', 'G', 'C', 'A'};

	const char* deltaA = deltas;
	const char* deltaC = deltas + 6;
	const char* deltaG = deltas + 12;
	const char* deltaT = deltas + 18;

	if (useConstDelta_)
	{
		rec_.sequence--;
		rec_.quality--;

		++rec_.sequenceLen;
		++rec_.qualityLen;
	}

	uchar symbol = seqStart_;
	rec_.sequence[0] = seqStart_;
	rec_.quality[0] = quaStart_;

	const char* dSelect = deltaA;

	for (uint32 k = 1; k < rec_.sequenceLen; ++k)
	{
		switch (symbol)
		{
			case 'A':	dSelect = deltaA;	break;
			case 'C':	dSelect = deltaC;	break;
			case 'G':	dSelect = deltaG;	break;
			case 'T':	dSelect = deltaT;	break;
			case 'N':	// use previous symbol matrix
			default:	break;
		}
		symbol = rec_.sequence[k];
		rec_.sequence[k] = (uchar) (std::find(dSelect, dSelect+6, rec_.sequence[k]) - dSelect) + '.';
	}
}


void IRecordsProcessor::InitializeStats()
{
	dnaStats.Clear();
	qualityStats.Clear();
	csStats.Clear();
	// need to clear ColorSpace stats?
}

void IRecordsProcessor::FinalizeStats()
{
	// process DNA stats
	dnaStats.symbolCount = 0;
	for (uint32 i = 0; i < DnaStats::MaxSymbolCount; ++i)
	{
		if (dnaStats.symbolFreqs[i] > 0)
		{
			dnaStats.symbols[i] = dnaStats.symbolCount++;
		}
	}

	// process quality stats
	qualityStats.symbolCount = 0;
	for (uint32 i = 0; i < QualityStats::MaxSymbolCount; ++i)
	{
		if (qualityStats.symbolFreqs[i] > 0)
		{
			qualityStats.symbols[i] = qualityStats.symbolCount++;
		}
	}
}

fq::FastqChecksum IRecordsProcessor::ProcessForward(FastqRecord *records_, uint64 recordsCount_, uint32 flags_)
{
	if (flags_ == fq::FastqChecksum::CALC_NONE)
	{
		for (uint64 i = 0; i < recordsCount_; ++i)
			ProcessForward(records_[i]);
		return fq::FastqChecksum();
	}

	fastqHasher.Reset();
	fastqHasher.SetFlags(flags_);
	for (uint64 i = 0; i < recordsCount_; ++i)
	{
		fq::FastqRecord& r = records_[i];
		fastqHasher.Update(r);
		ProcessForward(r);
	}

	return fastqHasher.GetChecksum();
}

fq::FastqChecksum IRecordsProcessor::ProcessBackward(FastqRecord *records_, uint64 recordsCount_, uint32 flags_)
{
	if (flags_ == fq::FastqChecksum::CALC_NONE)
	{
		for (uint64 i = 0; i < recordsCount_; ++i)
			ProcessBackward(records_[i]);
		return fq::FastqChecksum();
	}

	fastqHasher.Reset();
	fastqHasher.SetFlags(flags_);
	for (uint64 i = 0; i < recordsCount_; ++i)
	{
		fq::FastqRecord& r = records_[i];
		ProcessBackward(r);
		fastqHasher.Update(r);
	}

	return fastqHasher.GetChecksum();
}


// Lossless processor
//
LosslessRecordsProcessor::LosslessRecordsProcessor(uint32 qualityOffset_, bool colorSpace_)
	:	IRecordsProcessor(qualityOffset_, colorSpace_)
{
	std::fill(dnaToIndexTable, dnaToIndexTable + 128, +InvalidValue);
	std::fill(dnaFromIndexTable, dnaFromIndexTable + DnaStats::MaxSymbolCount, +InvalidValue);

	// normal symbols
	dnaToIndexTable[(int32)'A'] = 0;		dnaFromIndexTable[0]  = 'A';
	dnaToIndexTable[(int32)'G'] = 1;		dnaFromIndexTable[1]  = 'G';
	dnaToIndexTable[(int32)'C'] = 2;		dnaFromIndexTable[2]  = 'C';
	dnaToIndexTable[(int32)'T'] = 3;		dnaFromIndexTable[3]  = 'T';
	// amb codes
	dnaToIndexTable[(int32)'N'] = 4;		dnaFromIndexTable[4]  = 'N';
	dnaToIndexTable[(int32)'R'] = 5;		dnaFromIndexTable[5]  = 'R';
	dnaToIndexTable[(int32)'W'] = 6;		dnaFromIndexTable[6]  = 'W';
	dnaToIndexTable[(int32)'S'] = 7;		dnaFromIndexTable[7]  = 'S';
	dnaToIndexTable[(int32)'K'] = 8;		dnaFromIndexTable[8]  = 'K';
	dnaToIndexTable[(int32)'M'] = 9;		dnaFromIndexTable[9]  = 'M';
	dnaToIndexTable[(int32)'D'] = 10;		dnaFromIndexTable[10] = 'D';
	dnaToIndexTable[(int32)'V'] = 11;		dnaFromIndexTable[11] = 'V';
	dnaToIndexTable[(int32)'H'] = 12;		dnaFromIndexTable[12] = 'H';
	dnaToIndexTable[(int32)'B'] = 13;		dnaFromIndexTable[13] = 'B';
	dnaToIndexTable[(int32)'Y'] = 14;		dnaFromIndexTable[14] = 'Y';
	dnaToIndexTable[(int32)'X'] = 15;		dnaFromIndexTable[15] = 'X';
	dnaToIndexTable[(int32)'U'] = 16;		dnaFromIndexTable[16] = 'U';
	dnaToIndexTable[(int32)'.'] = 17;		dnaFromIndexTable[17] = '.';
	dnaToIndexTable[(int32)'-'] = 18;		dnaFromIndexTable[18] = '-';
}

void LosslessRecordsProcessor::ProcessForward(FastqRecord &rec_)
{
	uint32 seqLen = 0;
	uchar prevQSymbol = 255;
	uint32 curQThLen = 0;

	if (colorSpace)
	{
		ProcessFromColorSpace(rec_);
	}

	for (uint32 i = 0; i < rec_.sequenceLen; ++i)
	{
		ASSERT(dnaToIndexTable[rec_.sequence[i]] != DnaStats::EmptySymbol);
//		ASSERT(rec_.quality[i] >= qualityOffset && rec_.quality[i] - qualityOffset < 45);

		// transfom from symbol space to index space
		rec_.sequence[i] = dnaToIndexTable[rec_.sequence[i]];
		rec_.quality[i] -= qualityOffset;

		// check if the symbol differs from AGCT and whether can be transfered to quality stream
		if (rec_.sequence[i] > 3 && rec_.quality[i] < 7)
		{
			rec_.quality[i] += (uchar)(128 + (((uint32)rec_.sequence[i] - 3 + 1) << 3) - 16);
		}
		else
		{
			rec_.sequence[seqLen++] = rec_.sequence[i];

			// update DNA stats
			dnaStats.symbolFreqs[rec_.sequence[i]]++;
		}

		// update quality stats
		qualityStats.symbolFreqs[rec_.quality[i]]++;

		if (rec_.quality[i] != prevQSymbol)
			qualityStats.rleLength++;

		if (rec_.quality[i] != HashSymbolNormal)
			curQThLen = i;

		prevQSymbol = rec_.quality[i];
	}

	rec_.sequenceLen = seqLen;
	rec_.truncatedLen = curQThLen;
	rec_.truncatedLen += (uint16)(rec_.qualityLen > 0);

	// finalize quality stats
	if (prevQSymbol == HashSymbolNormal && qualityStats.rleLength > 0)
		qualityStats.rleLength--;

	qualityStats.rawLength += rec_.qualityLen;
	qualityStats.thLength += curQThLen;

	qualityStats.minLength = MIN(qualityStats.minLength, rec_.qualityLen);
	qualityStats.maxLength = MAX(qualityStats.maxLength, rec_.qualityLen);
}

void LosslessRecordsProcessor::ProcessBackward(FastqRecord &rec_)
{
	int32 seqi = rec_.sequenceLen-1;

	for (int32 i = rec_.qualityLen - 1; i >= 0; --i)
	{
		uint32 qval = rec_.quality[i];
		uint32 seqval = 0;

		// was ambiguous DNA symbol tramsfered to quality stream?
		if (qval >= 128)
		{
			seqval = (qval - 128 + 16)/8 + 3 - 1;
			ASSERT(seqval <= 18);
			qval &= 7;
		}
		else
		{
			seqval = rec_.sequence[seqi--];
		}

		// transfer from index space to symbol space
		ASSERT(seqval < DnaStats::MaxSymbolCount);
		rec_.sequence[i] = dnaFromIndexTable[seqval];
		rec_.quality[i] = qualityOffset + qval;
	}
	rec_.sequenceLen = rec_.qualityLen;

	if (colorSpace)
	{
		uchar seq0, qua0;
		if (csStats.constBeginSym)
		{
			seq0 = csStats.seqBegin;
			qua0 = csStats.quaBegin;
		}
		else
		{
			seq0 = rec_.sequence[0];
			qua0 = rec_.quality[0];
		}
		seq0 = dnaFromIndexTable[seq0];
		qua0 += qualityOffset;

		ProcessToColorSpace(rec_, seq0, qua0);
	}
}


LossyRecordsProcessor::LossyRecordsProcessor(uint32 qualityOffset_, bool colorSpace_)
	:	LosslessRecordsProcessor(qualityOffset_, colorSpace_)
{
	std::fill(qualityToIndexTable, qualityToIndexTable + 64, +InvalidValue);
	std::fill(qualityFromIndexTable, qualityFromIndexTable + 8, +InvalidValue);

	const uint32 ranges[] = {0, 2, 10, 20, 25, 30, 35, 40, 64};
	const uint32 qValues[] = {0, 6, 15, 22, 27, 33, 37, 40};
	const uint32 iValues[] = {0, 1, 2, 3, 4, 5, 6, 7};

	for (uint32 i = 0; i < 8; ++i)
	{
		uint32 r0 = ranges[i];
		uint32 r1 = ranges[i + 1];
		for (uint32 j = r0; j < r1; ++j)
		{
			qualityToIndexTable[j] = iValues[i];
		}
	}

	for (uint32 i = 0; i < 8; ++i)
	{
		qualityFromIndexTable[i] = qValues[i];
	}
}

void LossyRecordsProcessor::ProcessForward(FastqRecord &rec_)
{
	uchar prevQSymbol = 255;
	uint32 curQThLen = 0;

	if (colorSpace)
	{
		ProcessFromColorSpace(rec_);
	}

	uint32 seqLen = 0;
	for (uint32 i = 0; i < rec_.sequenceLen; ++i)
	{
		ASSERT(dnaToIndexTable[rec_.sequence[i]] != 255);
		ASSERT(rec_.quality[i] >= qualityOffset && rec_.quality[i] - qualityOffset < 45);

		// transfom from symbol to index space
		rec_.sequence[i] = dnaToIndexTable[rec_.sequence[i]];
		rec_.quality[i] = qualityToIndexTable[rec_.quality[i] - qualityOffset];

		// trim AMB code to N?
		if (rec_.sequence[i] >= 4)
		{
			// in most cases quality should be 0, if not - downgrade quality
			if (rec_.quality[i] != 0)
				rec_.quality[i] = 0;
		}
		else
		{
			// in most cases quality should be different from 0
			if (rec_.quality[i] == 0)
				rec_.quality[i] = 1;

			rec_.sequence[seqLen++] = rec_.sequence[i];

			dnaStats.symbolFreqs[rec_.sequence[i]]++;
		}

		qualityStats.symbolFreqs[rec_.quality[i]]++;

		if (rec_.quality[i] != prevQSymbol)
			qualityStats.rleLength++;

		if (rec_.quality[i] != HashSymbolNormal)
			curQThLen = i;

		prevQSymbol = rec_.quality[i];
	}

	rec_.sequenceLen = seqLen;

	rec_.sequenceLen = seqLen;
	rec_.truncatedLen = curQThLen;
	rec_.truncatedLen += (uint16)(rec_.qualityLen > 0);

	// finalize quality stats
	if (prevQSymbol == HashSymbolNormal && qualityStats.rleLength > 0)
		qualityStats.rleLength--;

	qualityStats.rawLength += rec_.qualityLen;
	qualityStats.thLength += curQThLen;

	qualityStats.minLength = MIN(qualityStats.minLength, rec_.qualityLen);
	qualityStats.maxLength = MAX(qualityStats.maxLength, rec_.qualityLen);
}

void LossyRecordsProcessor::ProcessBackward(FastqRecord &rec_)
{
	int32 seqi = rec_.sequenceLen-1;

	for (int32 i = rec_.qualityLen - 1; i >= 0; --i)
	{
		uint32 qval = rec_.quality[i];
		uint32 seqval = 0;

		// untransfer N ?
		if (qval == 0)
		{
			seqval = 4;
		}
		else
		{
			seqval = rec_.sequence[seqi--];
		}

		// from index to symbol
		rec_.sequence[i] = dnaFromIndexTable[seqval];
		rec_.quality[i] = qualityOffset + qualityFromIndexTable[qval];
	}
	rec_.sequenceLen = rec_.qualityLen;


	if (colorSpace)
	{
		uchar seq0, qua0;
		if (csStats.constBeginSym)
		{
			seq0 = csStats.seqBegin;
			qua0 = csStats.quaBegin;
		}
		else
		{
			seq0 = rec_.sequence[0];
			qua0 = rec_.quality[0];
		}
		seq0 = dnaFromIndexTable[seq0];
		qua0 += qualityOffset;

		ProcessToColorSpace(rec_, seq0, qua0);
	}
}

} // namespace comp

} // namespace dsrc
