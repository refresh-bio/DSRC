/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_BLOCKCOMPRESSOR
#define H_BLOCKCOMPRESSOR

#include "../include/dsrc/Globals.h"

#include <vector>

#include "Common.h"
#include "Fastq.h"
#include "RecordsProcessor.h"
#include "DnaModelerProxy.h"
#include "QualityModelerProxy.h"
#include "TagModeler.h"
#include "Crc32.h"

#include "huffman.h"

namespace dsrc
{

namespace comp
{

struct ChunkHeader
{
	uint64 recordsCount;
	uint64 chunkSize;
	uint32 flags;

	uint16 minQuaLength;
	uint16 maxQuaLength;

	bool csConstBeginSym;
	uchar csSeqBegin;
	uchar csQuaBegin;

	fq::FastqChecksum checksum;
	uint32 checksumFlags;

	ChunkHeader()
		:	recordsCount(0)
		,	chunkSize(0)
		,	flags(0)
		,	minQuaLength(-1)
		,	maxQuaLength(0)
		,	csConstBeginSym(false)
		,	csSeqBegin(0)
		,	csQuaBegin(0)
		,	checksumFlags(fq::FastqChecksum::CALC_NONE)
	{}
};


class BlockCompressor
{
public:
	BlockCompressor(const fq::FastqDatasetType& type_, const CompressionSettings& settings_);
	virtual ~BlockCompressor();
	
	void Store(core::BitMemoryWriter &memory_, fq::StreamsInfo& rawStreamsInfo_, fq::StreamsInfo& compStreamsInfo_, const fq::FastqDataChunk& chunk_);
	void Read(core::BitMemoryReader &memory_, fq::FastqDataChunk& chunk_);
	bool VerifyChecksum(core::BitMemoryReader &memory_, fq::FastqDataChunk& chunk_);

	void Reset();

protected:
	enum FastqBlockFlags
	{
		FLAG_DELTA_CONSTANT			= BIT(0),
		FLAG_VARIABLE_LENGTH		= BIT(1),
		FLAG_MIXED_FIELD_FORMATTING	= BIT(2)		// this should be handled by TagModelerProxy*
	};

	const fq::FastqDatasetType datasetType;
	const CompressionSettings compSettings;

	std::vector<fq::FastqRecord> records;

	ChunkHeader chunkHeader;

	IRecordsProcessor* recordsProcessor;
	TagModeler tagModeler;
	IDnaModelerProxy* dnaModeler;
	IQualityModeler* qualityModeler;

	void ParseRecords(const fq::FastqDataChunk& chunk_, fq::StreamsInfo& streamSizes_);

	void PreprocessRecords(uint32 checksumFlags_ = fq::FastqChecksum::CALC_NONE);
	void PostprocessRecords(uint32 checksumFlags_ = fq::FastqChecksum::CALC_NONE);

	void AnalyzeRecords();
	void AnalyzeMetaData(const DnaStats& dnaStats_, const QualityStats& qStats_, const ColorSpaceStats& csStats_);

	void StoreRecords(core::BitMemoryWriter &memory_, fq::StreamsInfo& streamInfo_);
	void ReadRecords(core::BitMemoryReader &memory_, fq::FastqDataChunk& chunk_);

	void StoreMetaData(core::BitMemoryWriter &memory_);
	void ReadMetaData(core::BitMemoryReader &memory_);

	void AnalyzeTags();
	void StoreTags(core::BitMemoryWriter &memory_);
	void ReadTags(core::BitMemoryReader &memory_, fq::FastqDataChunk& chunk_);

	void StoreDNA(core::BitMemoryWriter &memory_);
	void StoreQuality(core::BitMemoryWriter &memory_);

	void ReadDNA(core::BitMemoryReader &memory_);
	void ReadQuality(core::BitMemoryReader &memory_);
};

} // namespace comp

} // namespace dsrc

#endif

