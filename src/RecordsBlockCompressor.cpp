/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#include "RecordsBlockCompressor.h"
#include "FastqParser.h"

namespace dsrc
{

namespace ext
{

using namespace core;


void RecordsBlockCompressor::WriteNextRecord(const FastqRecord& rec_)
{
	std::vector<fq::FastqRecord>& records = compressor.records;

	if (recordsIdx > records.size())
		records.resize(records.size() + FastqRecordsIncr);

	// swap records data
	InsertNewRecord(rec_);
	recordsIdx++;
	recordsCount++;
}

void RecordsBlockCompressor::Flush(BitMemoryWriter &memory_)
{
	compressor.chunkHeader.recordsCount = recordsIdx;
	compressor.chunkHeader.rawChunkSize = rawChunkSize - 1;		// the last newline symbol in the block

	// from store()
	compressor.PreprocessRecords();

	compressor.AnalyzeRecords();

	fq::StreamsInfo info;
	compressor.StoreRecords(memory_, info);

	compressor.Reset();
	//


	fastqBuffer.size = 0;
	recordsIdx = 0;
	recordsCount = 0;
	rawChunkSize = 0;
}

bool RecordsBlockCompressor::ReadNextRecord(FastqRecord& rec_)
{
	if (recordsIdx >= recordsCount)
		return false;

	ExtractNextRecord(rec_);
	recordsIdx++;
	return true;
}

void RecordsBlockCompressor::Feed(BitMemoryReader &memory_)
{
	compressor.ReadRecords(memory_, fastqBuffer);
	ASSERT(fastqBuffer.size > 0);

	compressor.PostprocessRecords();
	recordsCount = compressor.chunkHeader.recordsCount;
	rawChunkSize = compressor.chunkHeader.rawChunkSize;

	recordsIdx = 0;
}

bool RecordsBlockCompressor::AnalyzeRecords(bool estimateQualityOffset_, bool &isColorSpace_, uint32 &qualityOffset_)
{
	return fq::FastqParser::Analyze(compressor.GetRecords(),
									recordsCount,
									estimateQualityOffset_,
									isColorSpace_,
									qualityOffset_);
}

void RecordsBlockCompressor::InsertNewRecord(const FastqRecord& rec_)
{
	uint64 rsize = RecordSize(rec_);
	std::vector<fq::FastqRecord>& records = compressor.records;

	if (fastqBuffer.size + rsize > fastqBuffer.data.Size())
	{
		fastqBuffer.data.Extend(fastqBuffer.data.Size() + (fastqBuffer.data.Size() / 2) + FastqBufferPadding, true);

		// rebind records to the new buffer
		byte* p = fastqBuffer.data.Pointer();
		for (uint64 i = 0; i < recordsIdx; ++i)
		{
			fq::FastqRecord& r = records[i];

			r.title = p;
			p += r.titleLen;
			r.sequence = p;
			p += r.sequenceLen;
			r.quality = p;
			p += r.qualityLen;
		}
	}

	if (recordsIdx + 1 > records.size())
		records.resize(records.size() + FastqRecordsIncr);

	fq::FastqRecord& r = records[recordsIdx];
	byte* p = fastqBuffer.data.Pointer() + fastqBuffer.size;

	// as this functionality comes directly from FastqParser,
	// here we need to apply TAG filter manually...
	const uint64 flags = compressor.GetCompressionSettings().tagPreserveFlags;
	if (flags != 0)
	{
		const char *fieldSeparators = " ._,=:/-#"; //9
		uint32 fieldNo = 0;
		uint32 fieldBeginPos = 0;
		uint32 bufferPos = 0;

		for (uint32 i = 0; i <= rec_.tag.length(); ++i)
		{
			if (!std::count(fieldSeparators, fieldSeparators + 10, rec_.tag[i]) && (i != rec_.tag.length()))
				continue;

			fieldNo++;

			if (BIT_ISSET(flags, fieldNo))
			{
				std::copy(rec_.tag.c_str() + fieldBeginPos, rec_.tag.c_str() + i + 1, p + bufferPos);
				bufferPos += (i + 1 - fieldBeginPos);
			}
			fieldBeginPos = i + 1;
		}

		// skip copying the separator after the last token
		if (bufferPos > 0 && bufferPos != fieldBeginPos)
			bufferPos -= 1;

		r.title = p;
		r.titleLen = bufferPos;
		p += bufferPos;
		fastqBuffer.size += bufferPos;

		// remember to update the raw FASTQ chunk size for decompression
		rsize -= rec_.tag.length() - bufferPos;
	}
	else
	{
		std::copy(rec_.tag.data(), rec_.tag.data() + rec_.tag.length(), p);
		r.title = p;
		r.titleLen = rec_.tag.length();
		p += rec_.tag.length();
		fastqBuffer.size += rec_.tag.length();
	}

	std::copy(rec_.sequence.data(), rec_.sequence.data() + rec_.sequence.length(), p);
	r.sequence = p;
	r.sequenceLen = rec_.sequence.length();
	p += rec_.sequence.length();
	fastqBuffer.size += rec_.sequence.length();

	std::copy(rec_.quality.data(), rec_.quality.data() + rec_.quality.length(), p);
	r.quality = p;
	r.qualityLen = rec_.quality.length();
	p += rec_.quality.length();
	fastqBuffer.size += rec_.quality.length();

	rawChunkSize += rsize;
}

uint64 RecordsBlockCompressor::RecordSize(const FastqRecord& r_)
{
	return r_.tag.length() + 1				// +1 as newline character
			+ r_.sequence.length() + 1
			+ r_.plus.length() + 1
			+ r_.quality.length() + 1;
}

void RecordsBlockCompressor::ExtractNextRecord(FastqRecord& rec_)
{
	fq::FastqRecord& r = compressor.records[recordsIdx];

	rec_.tag.assign(r.title, r.title + r.titleLen);
	rec_.sequence.assign(r.sequence, r.sequence + r.sequenceLen);
	if (compressor.datasetType.plusRepetition)
	{
		rec_.plus = rec_.tag;
		rec_.plus[0] = '+';
	}
	else if (rec_.plus.length() != 1)
	{
		rec_.plus.assign(1, '+');
		if (compressor.datasetType.plusRepetition)
		{
			rec_.plus.insert(rec_.plus.end(), rec_.tag.begin() + 1, rec_.tag.end());
		}
	}
	rec_.quality.assign(r.quality, r.quality + r.qualityLen);
}

} // namespace ext

} // namespace dsrc
