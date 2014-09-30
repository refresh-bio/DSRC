/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
*/

#include "BlockCompressorExt.h"

namespace dsrc
{

namespace wrap
{

using namespace core;

void BlockCompressorExt::WriteNextRecord(const FastqRecord& rec_)
{
	if (recordsIdx > records.size())
		records.resize(records.size() + FastqRecordsIncr);

	// swap records data
	InsertNewRecord(rec_);
	recordsIdx++;
}

void BlockCompressorExt::Flush(BitMemoryWriter &memory_)
{
	chunkHeader.recordsCount = recordsIdx;

	// from store()
	PreprocessRecords();

	AnalyzeRecords();

	fq::StreamsInfo info;
	StoreRecords(memory_, info);

	Reset();
	//

	fastqChunk.size = 0;
	recordsIdx = 0;
}

bool BlockCompressorExt::ReadNextRecord(FastqRecord& rec_)
{
	if (recordsIdx >= chunkHeader.recordsCount)
		return false;

	ExtractNextRecord(rec_);
	recordsIdx++;
	return true;
}

void BlockCompressorExt::Feed(BitMemoryReader &memory_)
{
	ReadRecords(memory_, fastqChunk);

	PostprocessRecords();

	recordsIdx = 0;
}

void BlockCompressorExt::InsertNewRecord(const FastqRecord& rec_)
{
	uint64 rsize = RecordSize(rec_);
	if (fastqChunk.size + rsize > fastqChunk.data.Size())
	{
		fastqChunk.data.Extend(fastqChunk.data.Size() + (fastqChunk.data.Size() / 2) + FastqBufferPadding, true);

		// rebind records to the new buffer
		byte* p = fastqChunk.data.Pointer();
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
	byte* p = fastqChunk.data.Pointer() + fastqChunk.size;

	std::copy(rec_.tag.data(), rec_.tag.data() + rec_.tag.length(), p);
	r.title = p;
	r.titleLen = rec_.tag.length();
	p += rec_.tag.length();
	fastqChunk.size += rec_.tag.length();

	std::copy(rec_.sequence.data(), rec_.sequence.data() + rec_.sequence.length(), p);
	r.sequence = p;
	r.sequenceLen = rec_.sequence.length();
	p += rec_.sequence.length();
	fastqChunk.size += rec_.sequence.length();

	std::copy(rec_.quality.data(), rec_.quality.data() + rec_.quality.length(), p);
	r.quality = p;
	r.qualityLen = rec_.quality.length();
	p += rec_.quality.length();
	fastqChunk.size += rec_.quality.length();

	chunkHeader.chunkSize += rsize;
}

uint64 BlockCompressorExt::RecordSize(const FastqRecord& r_)
{
	return r_.tag.length() + 1				// +1 as newline character
			+ r_.sequence.length() + 1
			+ r_.plus.length() + 1
			+ r_.quality.length() + 1;
}

void BlockCompressorExt::ExtractNextRecord(FastqRecord& rec_)
{
	fq::FastqRecord& r = records[recordsIdx];

	rec_.tag.assign(r.title, r.title + r.titleLen);
	rec_.sequence.assign(r.sequence, r.sequence + r.sequenceLen);
	if (datasetType.plusRepetition)
	{
		rec_.plus = rec_.tag;
		rec_.plus[0] = '+';
	}
	else if (rec_.plus.length() != 1)
	{
		rec_.plus.assign(1, '+');
		if (datasetType.plusRepetition)
		{
			rec_.plus.insert(rec_.plus.end(), rec_.tag.begin() + 1, rec_.tag.end());
		}
	}
	rec_.quality.assign(r.quality, r.quality + r.qualityLen);
}


} // namespace wrap

} // namespace dsrc
