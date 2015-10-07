/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef BLOCKCOMPRESSOREXT_H
#define BLOCKCOMPRESSOREXT_H

#include "../include/dsrc/Globals.h"
#include "../include/dsrc/FastqRecord.h"

#include "BlockCompressor.h"

namespace dsrc
{

namespace ext
{

#if 0

class BlockCompressorExt : public comp::BlockCompressor
{
public:
	BlockCompressorExt(const FastqDatasetType& type_, const comp::CompressionSettings& settings_)
		:	comp::BlockCompressor(type_, settings_)
		,	recordsIdx(0)
		,	fastqChunk(NULL)
	{
		//fastqChunk.data.Extend(DefaultFastqBufferSize + FastqBufferPadding);
	}

	void BindFastqChunkBuffer(fq::FastqDataChunk* chunk_)
	{
		ASSERT(chunk_ != NULL);
		fastqChunk = chunk_;
		fastqChunk->Reset();
	}

	uint64 ChunkSize()
	{
		return fastqChunk->size;
	}

	uint64 ChunkRecordsCount()
	{
		return chunkHeader.recordsCount;		// will be 0 when compressing
	}

	uint64 ChunkRecordsIdx() const
	{
		return recordsIdx;
	}

	void WriteNextRecord(const FastqRecord& rec_);
	void Flush(core::BitMemoryWriter &memory_);

	bool ReadNextRecord(FastqRecord& rec_);
	void Feed(core::BitMemoryReader &memory_);

private:
	static const uint64 FastqBufferPadding = 1 << 10 << 8;
	static const uint64 FastqRecordsIncr = 1 << 12;

	uint64 recordsIdx;
	fq::FastqDataChunk* fastqChunk;

	void InsertNewRecord(const FastqRecord& rec_);
	void ExtractNextRecord(FastqRecord& rec_);

	uint64 RecordSize(const FastqRecord& r_);
};

#endif


class RecordsBlockCompressor
{
public:
	RecordsBlockCompressor(comp::BlockCompressor& compressor_,
						   fq::FastqDataChunk& fastqBuffer_)
		:	compressor(compressor_)
		,	fastqBuffer(fastqBuffer_)
		,	recordsIdx(0)
		,	recordsCount(0)
		,	rawChunkSize(0)
	{}

	uint64 RawChunkSize()
	{
		return rawChunkSize;
	}

	uint64 ChunkRecordsCount()
	{
		return recordsCount;			// constant when decompressing
	}

	uint64 ChunkRecordsIdx() const
	{
		return recordsIdx;
	}

	void Reset()
	{
		recordsIdx = 0;
		recordsCount = 0;
		rawChunkSize = 0;
	}

	void WriteNextRecord(const FastqRecord& rec_);
	void Flush(core::BitMemoryWriter &memory_);

	bool ReadNextRecord(FastqRecord& rec_);
	void Feed(core::BitMemoryReader &memory_);

	bool AnalyzeRecords(bool estimateQualityOffset_, bool& isColorSpace_, uint32& qualityOffset_);

private:
	static const uint64 FastqBufferPadding = 1 << 10 << 8;
	static const uint64 FastqRecordsIncr = 1 << 12;

	comp::BlockCompressor& compressor;
	fq::FastqDataChunk& fastqBuffer;

	uint64 recordsIdx;
	uint64 recordsCount;
	uint64 rawChunkSize;

	void InsertNewRecord(const FastqRecord& rec_);
	void ExtractNextRecord(FastqRecord& rec_);

	uint64 RecordSize(const FastqRecord& r_);
};





} // namespace ext

} // namespace dsrc


#endif // BLOCKCOMPRESSOR_H
