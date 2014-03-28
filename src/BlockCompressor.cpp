/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "BlockCompressor.h"

#include <algorithm>
#include <cstring>
#include <ctime>

#include "BlockCompressor.h"
#include "BitMemory.h"
#include "FastqParser.h"

#include "utils.h"

namespace dsrc
{

namespace comp
{

using namespace core;
using namespace fq;

#if !defined(DEBUG) || (defined(DEBUG) && !DEBUG)

inline void CONTROL_CHECK_W(BitMemoryWriter& ) {};
inline void CONTROL_CHECK_R(BitMemoryReader& ) {};

#else

inline void CONTROL_CHECK_W(BitMemoryWriter& w_)
{
	uint32 pos = w_.Position();
	w_.PutWord(pos);
}

inline void CONTROL_CHECK_R(BitMemoryReader& r_)
{
	uint32 buf = r_.GetWord();
	ASSERT(r_.Position() == buf + 4);
}

#endif


// ********************************************************************************************
BlockCompressor::BlockCompressor(const FastqDatasetType& type_, const CompressionSettings& settings_)
	:	datasetType(type_)
	,	compSettings(settings_)
	,	recordsProcessor(NULL)
	,	dnaModeler(NULL)
	,	qualityModeler(NULL)
{
	records.resize(8 * 1024);

	if (settings_.lossy)
		recordsProcessor = new LossyRecordsProcessor(type_.qualityOffset, type_.colorSpace);
	else
		recordsProcessor = new LosslessRecordsProcessor(type_.qualityOffset, type_.colorSpace);

	if (settings_.dnaOrder == 0)
		dnaModeler = new DnaNormalModelerProxy();
	else
		dnaModeler = new DnaOrderModelerProxy(settings_.dnaOrder);

	if (settings_.qualityOrder > 0)
	{
		if (settings_.lossy)
			qualityModeler = new QualityOrderModelerProxyLossy(settings_.qualityOrder);
		else
			qualityModeler = new QualityOrderModelerProxyLossless(settings_.qualityOrder);
	}
	else
	{
		qualityModeler = new QualityNormalModelerProxy(settings_.lossy);
	}

	if (settings_.calculateCrc32)
	{
		if (settings_.tagPreserveFlags == CompressionSettings::DefaultTagPreserveFlags)
			chunkHeader.checksumFlags |= fq::FastqChecksum::CALC_TAG;

		chunkHeader.checksumFlags |= fq::FastqChecksum::CALC_SEQUENCE;

		if (!settings_.lossy)
			chunkHeader.checksumFlags |= fq::FastqChecksum::CALC_QUALITY;
	}
}

// ********************************************************************************************
BlockCompressor::~BlockCompressor()
{
	delete qualityModeler;
	delete dnaModeler;
	delete recordsProcessor;
}

// ********************************************************************************************
void BlockCompressor::Reset()
{
	chunkHeader.flags = 0;
	chunkHeader.recordsCount = 0;
}

void BlockCompressor::ParseRecords(const FastqDataChunk& chunk_)
{
	if (compSettings.tagPreserveFlags != 0)
	{
		FastqParserExt parser;
		uint64 size = parser.ParseFrom(chunk_, records, chunkHeader.recordsCount, compSettings.tagPreserveFlags);
		ASSERT(size <= chunk_.size);
		chunkHeader.chunkSize = size;
	}
	else
	{
		FastqParser parser;
		uint64 size = parser.ParseFrom(chunk_, records, chunkHeader.recordsCount);
		ASSERT(size <= chunk_.size);
		chunkHeader.chunkSize = size;
	}

	ASSERT(chunkHeader.recordsCount > 0);
	ASSERT(records.size() >= chunkHeader.recordsCount);
}


void BlockCompressor::PreprocessRecords(uint32 checksumFlags_)
{
	recordsProcessor->InitializeStats();
	fq::FastqChecksum checksum = recordsProcessor->ProcessForward(records.data(), chunkHeader.recordsCount, checksumFlags_);
	recordsProcessor->FinalizeStats();

	if (checksumFlags_ != fq::FastqChecksum::CALC_NONE)
	{
		chunkHeader.checksum = checksum;
	}
}

void BlockCompressor::PostprocessRecords(uint32 checksumFlags_)
{
	if (datasetType.colorSpace)
	{
		ColorSpaceStats stats;
		stats.constBeginSym = chunkHeader.csConstBeginSym;
		stats.seqBegin = chunkHeader.csSeqBegin;
		stats.quaBegin = chunkHeader.csQuaBegin;
		recordsProcessor->SetColorSpaceStats(stats);
	}
	fq::FastqChecksum checksum = recordsProcessor->ProcessBackward(records.data(), chunkHeader.recordsCount, checksumFlags_);

	if (checksumFlags_ != fq::FastqChecksum::CALC_NONE)
	{
		chunkHeader.checksum = checksum;
	}
}

void BlockCompressor::AnalyzeRecords()
{	
	AnalyzeMetaData(recordsProcessor->GetDnaStats(), recordsProcessor->GetQualityStats(), recordsProcessor->GetColorSpaceStats());

	AnalyzeTags();

	dnaModeler->ProcessStats(recordsProcessor->GetDnaStats());

	qualityModeler->ProcessStats(recordsProcessor->GetQualityStats());
}


void BlockCompressor::AnalyzeMetaData(const DnaStats& , const QualityStats& qStats_, const ColorSpaceStats& csStats_)
{
	chunkHeader.maxQuaLength = qStats_.maxLength;
	chunkHeader.minQuaLength = qStats_.minLength;
	chunkHeader.csConstBeginSym = csStats_.constBeginSym;

	if (datasetType.colorSpace && csStats_.constBeginSym)
	{
		chunkHeader.flags |= FLAG_DELTA_CONSTANT;

		chunkHeader.csSeqBegin = records[0].sequence[0];
		chunkHeader.csQuaBegin = records[0].quality[0];

		chunkHeader.maxQuaLength -= 1;
		chunkHeader.minQuaLength -= 1;
	}

	if (chunkHeader.maxQuaLength != chunkHeader.minQuaLength)
	{
		chunkHeader.flags |= FLAG_VARIABLE_LENGTH;
	}
}

// ********************************************************************************************
void BlockCompressor::Store(BitMemoryWriter &memory_, const FastqDataChunk &chunk_)
{
	ParseRecords(chunk_);

	PreprocessRecords(chunkHeader.checksumFlags);

	AnalyzeRecords();

	StoreRecords(memory_);

	Reset();
}
void BlockCompressor::StoreRecords(BitMemoryWriter &memory_)
{
	CONTROL_CHECK_W(memory_);
	StoreMetaData(memory_);

	CONTROL_CHECK_W(memory_);
	StoreTags(memory_);

	CONTROL_CHECK_W(memory_);
	StoreQuality(memory_);

	CONTROL_CHECK_W(memory_);
	StoreDNA(memory_);

	CONTROL_CHECK_W(memory_);
}

// ********************************************************************************************
void BlockCompressor::Read(BitMemoryReader &memory_, FastqDataChunk &chunk)
{
	ReadRecords(memory_, chunk);

	PostprocessRecords(fq::FastqChecksum::CALC_NONE);

	Reset();
}

void BlockCompressor::ReadRecords(BitMemoryReader &memory_, FastqDataChunk &chunk_)
{
	CONTROL_CHECK_R(memory_);
	ReadMetaData(memory_);

	// extend chunk if necessary
	//
	chunkHeader.chunkSize += 1; // +1 for the last '\n'
	chunk_.size = chunkHeader.chunkSize;

	if (chunk_.data.Size() < chunkHeader.chunkSize)
	{
		chunk_.data.Extend(chunkHeader.chunkSize + MEM_EXTENSION_FACTOR(chunkHeader.chunkSize));
	}

	CONTROL_CHECK_R(memory_);
	ReadTags(memory_, chunk_);

	CONTROL_CHECK_R(memory_);
	ReadQuality(memory_);

	CONTROL_CHECK_R(memory_);
	ReadDNA(memory_);

	CONTROL_CHECK_R(memory_);
}

void BlockCompressor::ReadMetaData(BitMemoryReader &memory_)
{
	chunkHeader.recordsCount = memory_.GetWord();
	chunkHeader.maxQuaLength = memory_.GetWord();

	chunkHeader.flags = memory_.GetWord();
	ASSERT(chunkHeader.flags < 1 << 8);

	chunkHeader.chunkSize = memory_.GetWord();

	// setup records
	//
	if (records.size() < chunkHeader.recordsCount)
	{
		records.resize(chunkHeader.recordsCount + REC_EXTENSION_FACTOR(chunkHeader.recordsCount));
	}

	if ((chunkHeader.flags & FLAG_VARIABLE_LENGTH) != 0)
	{
		chunkHeader.minQuaLength = memory_.GetWord();
		ASSERT(chunkHeader.maxQuaLength >= chunkHeader.minQuaLength);
	}
	else
	{
		chunkHeader.minQuaLength = chunkHeader.maxQuaLength;
	}

	if (datasetType.colorSpace)
	{
		chunkHeader.csConstBeginSym = (chunkHeader.flags & FLAG_DELTA_CONSTANT) != 0;
		if (chunkHeader.csConstBeginSym)
		{
			chunkHeader.csSeqBegin = memory_.GetByte();
			chunkHeader.csQuaBegin = memory_.GetByte();
		}
	}

	if (compSettings.calculateCrc32)
	{
		if (compSettings.tagPreserveFlags == CompressionSettings::DefaultTagPreserveFlags)
		{
			chunkHeader.checksum.tag = memory_.GetWord();
			ASSERT(chunkHeader.checksum.tag != 0);
		}

		chunkHeader.checksum.sequence = memory_.GetWord();
		ASSERT(chunkHeader.checksum.sequence != 0);

		if (!compSettings.lossy)
		{
			chunkHeader.checksum.quality = memory_.GetWord();
			ASSERT(chunkHeader.checksum.quality != 0);
		}
	}

	memory_.FlushInputWordBuffer();
}


void BlockCompressor::AnalyzeTags()
{
	bool cs_reduce_lens = datasetType.colorSpace && chunkHeader.csConstBeginSym;

	tagModeler.InitializeFieldsStats(records[0]);

	for (uint32 j = 0; j < chunkHeader.recordsCount; ++j)
	{
		FastqRecord &rec = records[j];

		tagModeler.UpdateFieldsStats(rec);

		//
		// this should be logically split
		//

		// 2nd pass:
		//
		if (cs_reduce_lens)
		{
			ASSERT(rec.qualityLen > 1);
			ASSERT(rec.sequenceLen > 1);

			rec.sequence++;
			rec.quality++;

			rec.qualityLen -= 1;
			rec.sequenceLen -= 1;

			//ASSERT(rec.truncatedLen > 0);		// this can be buggy in case '#######...'
			if (rec.truncatedLen > 0)
				rec.truncatedLen -= 1;
		}
	}

	tagModeler.FinalizeFieldsStats();
}


void BlockCompressor::StoreMetaData(BitMemoryWriter &memory_)
{
	memory_.PutWord(chunkHeader.recordsCount);
	memory_.PutWord(chunkHeader.maxQuaLength);
	memory_.PutWord(chunkHeader.flags);
	memory_.PutWord(chunkHeader.chunkSize);

	if ((chunkHeader.flags & FLAG_VARIABLE_LENGTH) != 0)
	{
		memory_.PutWord(chunkHeader.minQuaLength);
	}

	if (datasetType.colorSpace)
	{
		if ((chunkHeader.flags & FLAG_DELTA_CONSTANT) != 0)
		{
			memory_.PutByte(chunkHeader.csSeqBegin);
			memory_.PutByte(chunkHeader.csQuaBegin);
		}
	}

	if (compSettings.calculateCrc32)
	{
		if (compSettings.tagPreserveFlags == CompressionSettings::DefaultTagPreserveFlags)
		{
			ASSERT(chunkHeader.checksum.tag != 0);
			memory_.PutWord(chunkHeader.checksum.tag);
		}

		ASSERT(chunkHeader.checksum.sequence != 0);
		memory_.PutWord(chunkHeader.checksum.sequence);

		if (!compSettings.lossy)
		{
			ASSERT(chunkHeader.checksum.quality != 0);
			memory_.PutWord(chunkHeader.checksum.quality);
		}
	}

	memory_.FlushPartialWordBuffer();
}

// ********************************************************************************************

void BlockCompressor::StoreDNA(BitMemoryWriter &memory_)
{
	dnaModeler->Encode(memory_, records.data(), chunkHeader.recordsCount);
}


void BlockCompressor::StoreQuality(BitMemoryWriter &memory_)
{
	qualityModeler->Encode(memory_, records.data(), chunkHeader.recordsCount);
}


// ********************************************************************************************
void BlockCompressor::StoreTags(BitMemoryWriter &memory_)
{
	const uint32 lenBits = core::bit_length(chunkHeader.maxQuaLength - chunkHeader.minQuaLength);
	const bool isVariableLen = lenBits > 0;

	tagModeler.StartEncoding(memory_);

	// store record title info + some meta-data
	//
	for (uint32 i = 0; i < chunkHeader.recordsCount; ++i)
	{
		FastqRecord &rec = records[i];

		tagModeler.EncodeNextFields(memory_, rec);

		// save other meta info
		//
		if (isVariableLen)
		{
			memory_.PutBits(rec.qualityLen - chunkHeader.minQuaLength, lenBits);
		}
	}

	tagModeler.FinishEncoding(memory_);

	//memory_.FlushPartialWordBuffer();
}


// ********************************************************************************************
void BlockCompressor::ReadDNA(BitMemoryReader &memory_)
{
	dnaModeler->Decode(memory_, records.data(), chunkHeader.recordsCount);
}

void BlockCompressor::ReadQuality(BitMemoryReader &memory_)
{
	qualityModeler->Decode(memory_, records.data(), chunkHeader.recordsCount);
}

// ********************************************************************************************
void BlockCompressor::ReadTags(BitMemoryReader &memory_, FastqDataChunk& fq_chunk_)
{
	uchar* chunkBegin = fq_chunk_.data.Pointer();
	uint32 bufPos = 0;

	memory_.FlushInputWordBuffer();

	const uint32 lenBits = core::bit_length(chunkHeader.maxQuaLength - chunkHeader.minQuaLength);
	const bool isVariableLen =  lenBits > 0;

	const bool cs_const_delta_enc = datasetType.colorSpace && (chunkHeader.flags & FLAG_DELTA_CONSTANT) != 0;

	tagModeler.StartDecoding(memory_);

	for (uint32 i = 0; i < chunkHeader.recordsCount; ++i)
	{
		ASSERT(bufPos < fq_chunk_.data.Size());

		FastqRecord& cur_rec = records[i];
		cur_rec.titleLen = 0;
		cur_rec.title = chunkBegin + bufPos;

		tagModeler.DecodeNextFields(memory_, cur_rec);

		// parsed the title, now bind the rest of pointers to buffer
		//
		bufPos += cur_rec.titleLen;		// title
		chunkBegin[bufPos++] = '\n';


		//
		// this can be logically split
		//


		// read here some record bit flags -- quality flags
		//
		if (isVariableLen)
		{
			cur_rec.qualityLen = memory_.GetBits(lenBits) + chunkHeader.minQuaLength;
		}
		else
		{
			cur_rec.qualityLen = chunkHeader.maxQuaLength;
		}

		cur_rec.sequenceLen = cur_rec.qualityLen;

		cur_rec.sequence = chunkBegin + bufPos;
		bufPos += cur_rec.sequenceLen;
		if (cs_const_delta_enc)
		{
			cur_rec.sequence++;
			bufPos++;
		}
		chunkBegin[bufPos++] = '\n';
		
		chunkBegin[bufPos++] = '+';
		if (datasetType.plusRepetition)
		{
			std::copy(cur_rec.title + 1, cur_rec.title + cur_rec.titleLen, chunkBegin + bufPos);
			bufPos += cur_rec.titleLen - 1;
		}
		chunkBegin[bufPos++] = '\n';

		cur_rec.quality = chunkBegin + bufPos;
		bufPos += cur_rec.qualityLen;
		if (cs_const_delta_enc)
		{
			cur_rec.quality++;
			bufPos++;
		}
		ASSERT(bufPos < fq_chunk_.size);
		chunkBegin[bufPos++] = '\n';
	}

	tagModeler.FinishDecoding(memory_);
}

bool BlockCompressor::VerifyChecksum(BitMemoryReader &memory_, FastqDataChunk &chunk)
{
	ASSERT(compSettings.calculateCrc32);

	ReadRecords(memory_, chunk);

	fq::FastqChecksum blockCrc32 = chunkHeader.checksum;
	PostprocessRecords(chunkHeader.checksumFlags);

	Reset();

	bool valid = true;
	if (compSettings.tagPreserveFlags == CompressionSettings::DefaultTagPreserveFlags)
		valid &= blockCrc32.tag == chunkHeader.checksum.tag;
	valid &= blockCrc32.sequence == chunkHeader.checksum.sequence;
	if (!compSettings.lossy)
		valid &= blockCrc32.quality == chunkHeader.checksum.quality;
	return valid;
}

} // namespace comp

} // namespace dsrc
