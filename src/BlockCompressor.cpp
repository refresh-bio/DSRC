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


BlockCompressor::~BlockCompressor()
{
	delete qualityModeler;
	delete dnaModeler;
	delete recordsProcessor;
}


void BlockCompressor::Reset()
{
	chunkHeader.flags = 0;
	chunkHeader.recordsCount = 0;
}


void BlockCompressor::ParseRecords(const FastqDataChunk& chunk_, StreamsInfo& streamInfo_)
{
	if (compSettings.tagPreserveFlags != 0)
	{
		FastqParserExt parser;
		uint64 size = parser.ParseFrom(chunk_, records, chunkHeader.recordsCount,
									   streamInfo_, compSettings.tagPreserveFlags);
		ASSERT(size <= chunk_.size);
		chunkHeader.chunkSize = size;
	}
	else
	{
		FastqParser parser;
		uint64 size = parser.ParseFrom(chunk_, records, chunkHeader.recordsCount, streamInfo_);
		ASSERT(size <= chunk_.size);
		chunkHeader.chunkSize = size;
	}

	ASSERT(chunkHeader.recordsCount > 0);
	ASSERT(records.size() >= chunkHeader.recordsCount);
	const uint64 rawStreamSize = streamInfo_.sizes[StreamsInfo::TagStream]
							   + streamInfo_.sizes[StreamsInfo::DnaStream]
							   + streamInfo_.sizes[StreamsInfo::QualityStream]
							   + chunkHeader.recordsCount * 5 - 1;		// add newlines and pluses
	ASSERT(chunkHeader.chunkSize >= rawStreamSize);
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


void BlockCompressor::Store(BitMemoryWriter &memory_, StreamsInfo& rawStreamInfo_, StreamsInfo& compStreamInfo_,
							const FastqDataChunk &chunk_)
{
	ParseRecords(chunk_, rawStreamInfo_);

	PreprocessRecords(chunkHeader.checksumFlags);

	AnalyzeRecords();

	StoreRecords(memory_, compStreamInfo_);

	Reset();
}


void BlockCompressor::StoreRecords(BitMemoryWriter &memory_, StreamsInfo& streamInfo_)
{
	uint64 pos = memory_.Position();

	// store meta data
	//
	CONTROL_CHECK_W(memory_);
	StoreMetaData(memory_);

	streamInfo_.sizes[StreamsInfo::MetaStream] = memory_.Position() - pos;
	pos = memory_.Position();

	// store tags
	//
	CONTROL_CHECK_W(memory_);
	StoreTags(memory_);

	streamInfo_.sizes[StreamsInfo::TagStream] = memory_.Position() - pos;
	pos = memory_.Position();

	// store quality
	//
	CONTROL_CHECK_W(memory_);
	StoreQuality(memory_);

	streamInfo_.sizes[StreamsInfo::QualityStream] = memory_.Position() - pos;
	pos = memory_.Position();

	// store dna
	//
	CONTROL_CHECK_W(memory_);
	StoreDNA(memory_);

	streamInfo_.sizes[StreamsInfo::DnaStream] = memory_.Position() - pos;

	CONTROL_CHECK_W(memory_);
}


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

	TagAnalyzer* analyzer = tagModeler.GetAnalyzer();

	analyzer->InitializeFieldsStats(records[0]);

	for (uint32 j = 0; j < chunkHeader.recordsCount; ++j)
	{
		FastqRecord &rec = records[j];

		analyzer->UpdateFieldsStats(rec);

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

	analyzer->FinalizeFieldsStats();

	if (analyzer->GetStats().mixedFormatting)
		chunkHeader.flags |= FLAG_MIXED_FIELD_FORMATTING;
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


void BlockCompressor::StoreDNA(BitMemoryWriter &memory_)
{
	dnaModeler->Encode(memory_, records.data(), chunkHeader.recordsCount);
}


void BlockCompressor::StoreQuality(BitMemoryWriter &memory_)
{
	qualityModeler->Encode(memory_, records.data(), chunkHeader.recordsCount);
}


void BlockCompressor::StoreTags(BitMemoryWriter &memory_)
{
	ITagEncoder* encoder = NULL;
	if ((chunkHeader.flags & FLAG_MIXED_FIELD_FORMATTING) != 0)
		encoder = tagModeler.SelectEncoder(TagModeler::TagRawHuffman);
	else
		encoder = tagModeler.SelectEncoder(TagModeler::TagTokenizeHuffman);

	const uint32 lenBits = core::bit_length(chunkHeader.maxQuaLength - chunkHeader.minQuaLength);
	const bool isVariableLen = lenBits > 0;

	encoder->StartEncoding(memory_, &tagModeler.GetAnalyzer()->GetStats());

	// store record title info + some meta-data
	//
	for (uint32 i = 0; i < chunkHeader.recordsCount; ++i)
	{
		const FastqRecord &rec = records[i];

		encoder->EncodeNextFields(memory_, rec);

		// save other meta info
		//
		if (isVariableLen)
		{
			memory_.PutBits(rec.qualityLen - chunkHeader.minQuaLength, lenBits);
		}
	}

	encoder->FinishEncoding(memory_);
}


void BlockCompressor::ReadDNA(BitMemoryReader &memory_)
{
	dnaModeler->Decode(memory_, records.data(), chunkHeader.recordsCount);
}


void BlockCompressor::ReadQuality(BitMemoryReader &memory_)
{
	qualityModeler->Decode(memory_, records.data(), chunkHeader.recordsCount);
}


void BlockCompressor::ReadTags(BitMemoryReader &memory_, FastqDataChunk& fqChunk_)
{
	ITagDecoder* decoder = NULL;

	if ((chunkHeader.flags & FLAG_MIXED_FIELD_FORMATTING) != 0)
		decoder = tagModeler.SelectDecoder(TagModeler::TagRawHuffman);
	else
		decoder = tagModeler.SelectDecoder(TagModeler::TagTokenizeHuffman);

	uchar* chunkBegin = fqChunk_.data.Pointer();
	uint32 bufPos = 0;

	const uint32 lenBits = core::bit_length(chunkHeader.maxQuaLength - chunkHeader.minQuaLength);
	const bool isVariableLen =  lenBits > 0;
	const bool csConstDeltaEncode = datasetType.colorSpace && (chunkHeader.flags & FLAG_DELTA_CONSTANT) != 0;

	decoder->StartDecoding(memory_);

	for (uint32 i = 0; i < chunkHeader.recordsCount; ++i)
	{
		ASSERT(bufPos < fqChunk_.data.Size());

		FastqRecord& curRec = records[i];
		curRec.titleLen = 0;
		curRec.title = chunkBegin + bufPos;

		decoder->DecodeNextFields(memory_, curRec);

		bufPos += curRec.titleLen;		// title
		chunkBegin[bufPos++] = '\n';

		// code below should be logically split, but to avoid another loop
		// throught the records we are doing it here

		if (isVariableLen)
			curRec.qualityLen = memory_.GetBits(lenBits) + chunkHeader.minQuaLength;
		else
			curRec.qualityLen = chunkHeader.maxQuaLength;

		curRec.sequenceLen = curRec.qualityLen;

		curRec.sequence = chunkBegin + bufPos;
		bufPos += curRec.sequenceLen;
		if (csConstDeltaEncode)
		{
			curRec.sequence++;
			bufPos++;
		}
		chunkBegin[bufPos++] = '\n';

		chunkBegin[bufPos++] = '+';
		if (datasetType.plusRepetition)
		{
			std::copy(curRec.title + 1, curRec.title + curRec.titleLen, chunkBegin + bufPos);
			bufPos += curRec.titleLen - 1;
		}
		chunkBegin[bufPos++] = '\n';

		curRec.quality = chunkBegin + bufPos;
		bufPos += curRec.qualityLen;
		if (csConstDeltaEncode)
		{
			curRec.quality++;
			bufPos++;
		}
		ASSERT(bufPos < fqChunk_.size);
		chunkBegin[bufPos++] = '\n';
	}

	decoder->FinishDecoding(memory_);
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
