/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include "FastqParser.h"

namespace dsrc
{

namespace fq
{

FastqParser::FastqParser()
	:	buffer(NULL)
	,	memory(NULL)
	,	memoryPos(0)
	,	memorySize(0)
	,	skippedBytes(0)
{}


bool FastqParser::Analyze(const FastqDataChunk& chunk_, FastqDatasetType& header_, bool estimateQualityOffset_)
{
	uchar minQuality = (uchar)-1;
	uchar maxQuality = 0;

	memory = (byte*)chunk_.data.Pointer();
	memoryPos = 0;
	memorySize = chunk_.size;

	header_.colorSpace = false;
	header_.plusRepetition = false;
	//header_.qualityOffset = 0;

	uint32 recCount = 0;
	while (memoryPos < memorySize)
	{
		// read stuff
		//
		byte* title = memory + memoryPos;
		uint32 titleLen = SkipLine();
		if (titleLen == 0 || title[0] != '@')
			break;

		byte* sequence = memory + memoryPos;
		uint32 seqLen = SkipLine();
		if (seqLen == 0)
			break;

		byte* plus = memory + memoryPos;
		bool plusRep = SkipLine() > 1;
		if (plus[0] != '+')
			break;

		if (estimateQualityOffset_)
		{
			byte* qua = memory + memoryPos;
			uint32 quaLen = SkipLine();

			for (uint32 i = 0; i < quaLen; ++i)
			{
				minQuality = MIN(minQuality, qua[i]);
				maxQuality = MAX(maxQuality, qua[i]);
			}
		}
		else
		{
			if (SkipLine() == 0)	// read quality
				break;
		}


		// analyze stuff
		//
		bool colorEnc = (sequence[1] >= '0' && sequence[1] <= '3') || sequence[1] == '.';
		if (recCount != 0)
		{
			if (header_.colorSpace != colorEnc)
			{
				//throw std::runtime_error("Inconsistent sequence (ColorSpace / Normal) format");
				return false;
			}

			if (header_.colorSpace)
			{
				//throw std::runtime_error("Invalid color-space format");
				if (sequence[0] >= '0' && sequence[0] <= '3')
					return false;
			}

			if (header_.plusRepetition != plusRep)
			{
				//throw std::runtime_error("Inconsistent plus field format");
				return false;
			}
		}
		else
		{
			header_.plusRepetition = plusRep;
			header_.colorSpace = colorEnc;
		}

		recCount++;
	}

	if (estimateQualityOffset_)
	{
		// standard quality scores
		if (maxQuality <= 74)
		{
			if (minQuality >= 33)
				header_.qualityOffset = 33;				// standard Sanger / Illumina 1.8+
		}
		else if (maxQuality <= 105)
		{
			if (minQuality >= 64)
				header_.qualityOffset = 64;				// Illumina 1.3-1.8
			else if (minQuality >= 59)
				header_.qualityOffset = 59;				// Solexa
		}

		// check non-standard
		if (header_.qualityOffset == 0)
		{
			if (minQuality >= 33)
				header_.qualityOffset = 33;
			else
				return false;
		}
	}

	return recCount > 1;
}

uint64 FastqParser::ParseFrom(const FastqDataChunk& chunk_, std::vector<FastqRecord>& records_, uint64& rec_count_, StreamsInfo& streamsInfo_)
{
	ASSERT(buffer == NULL);

	memory = (byte*)chunk_.data.Pointer();
	memoryPos = 0;
	memorySize = chunk_.size;

	streamsInfo_.Clear();
	rec_count_ = 0;
	while (memoryPos < memorySize && ReadNextRecord(records_[rec_count_]))
	{
		const FastqRecord& rec = records_[rec_count_];
		streamsInfo_.sizes[StreamsInfo::TagStream] += rec.titleLen;
		streamsInfo_.sizes[StreamsInfo::DnaStream] += rec.sequenceLen;
		streamsInfo_.sizes[StreamsInfo::QualityStream] += rec.qualityLen;

		rec_count_++;
		if (records_.size() < rec_count_ + 1)
			records_.resize(records_.size() + REC_EXTENSION_FACTOR(records_.size()));
	}
	ASSERT(rec_count_ > 0);

	return chunk_.size - skippedBytes;
}


uint64 FastqParserExt::ParseFrom(const FastqDataChunk &chunk_, std::vector<FastqRecord> &records_, uint64 &rec_count_,
								 StreamsInfo& streamsInfo_, uint64 tagPreserveFlags_)
{
	ASSERT(buffer == NULL);
	ASSERT(tagPreserveFlags_ != 0);

	uchar tagBuffer[MaxTagBufferSize];
	totalBytesCut = 0;

	memory = (byte*)chunk_.data.Pointer();
	memoryPos = 0;
	memorySize = chunk_.size;

	rec_count_ = 0;
	while (memoryPos < memorySize && ReadNextRecord(records_[rec_count_], tagBuffer, tagPreserveFlags_))
	{
		const FastqRecord& rec = records_[rec_count_];
		streamsInfo_.sizes[StreamsInfo::TagStream] += rec.titleLen;
		streamsInfo_.sizes[StreamsInfo::DnaStream] += rec.sequenceLen;
		streamsInfo_.sizes[StreamsInfo::QualityStream] += rec.qualityLen;

		rec_count_++;
		if (records_.size() < rec_count_ + 1)
			records_.resize(records_.size() + REC_EXTENSION_FACTOR(records_.size()));
	}
	ASSERT(rec_count_ > 0);
	ASSERT(chunk_.size >= totalBytesCut + skippedBytes);

	return chunk_.size - totalBytesCut - skippedBytes;
}

bool FastqParserExt::ReadNextRecord(FastqRecord &rec_, uchar *tagBuffer_, uint64 tagPreserveFlags_)
{
	const char *fieldSeparators = " ._,=:/-#"; //9
	//const std::vector<uchar> separators(fieldSeparators, fieldSeparators + 9 + 1);

	if (memoryPos == memorySize)
		return false;

	rec_.title = memory + memoryPos;
	rec_.titleLen = SkipLine();
	if (rec_.titleLen == 0 || rec_.title[0] != '@')
		return false;

	ASSERT(rec_.titleLen <= MaxTagBufferSize);

	uint32 fieldNo = 0;
	uint32 fieldBeginPos = 0;
	uint32 bufferPos = 0;
	for (uint32 i = 0; i <= rec_.titleLen; ++i)
	{
		if (!std::count(fieldSeparators, fieldSeparators + 10, rec_.title[i]) && (i != rec_.titleLen))
			continue;

		fieldNo++;

		if (BIT_ISSET(tagPreserveFlags_, fieldNo))
		{
			std::copy(rec_.title + fieldBeginPos, rec_.title + i + 1, tagBuffer_ + bufferPos);
			bufferPos += (i + 1 - fieldBeginPos);
		}
		fieldBeginPos = i + 1;
	}

	ASSERT(rec_.titleLen >= bufferPos);
	totalBytesCut += rec_.titleLen - bufferPos;

	if (bufferPos > 0)
	{
		std::copy(tagBuffer_, tagBuffer_ + bufferPos, rec_.title);
	}
	rec_.titleLen = bufferPos;


	rec_.sequence = memory + memoryPos;
	rec_.sequenceLen = SkipLine();

	// read plus
	uint32 plusLen = SkipLine();

	rec_.quality = memory + memoryPos;
	rec_.qualityLen = SkipLine();

	return (plusLen > 0 && rec_.sequenceLen == rec_.qualityLen);
}

} // namespace fq

} // namespace dsrc
