/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
#ifndef H_DSRCFILE
#define H_DSRCFILE

#include "../include/dsrc/Globals.h"

#include <vector>

#include "Common.h"
#include "FileStream.h"
#include "Fastq.h"

namespace dsrc
{

namespace comp
{

struct DsrcFileHeader
{
	static const uchar DummyByteValue			= 0xAA;
	static const uint32	ReservedBytes			= 8;
	static const uint32 HeaderSize				= 4 + ReservedBytes + 3*8 + 4;

	static const uint32 VersionMajor = 2;
	static const uint32 VersionMinor = 0;
	static const uint32 VersionRev = 2;

	uchar	dummyByte;
	uchar	versionMajor;
	uchar	versionMinor;
	uchar	versionRev;
	uint32	footerSize;

	uint64	footerOffset;
	uint64	recordsCount;
	uint64	blockCount;

	uchar	reserved[ReservedBytes];
};

struct DsrcFileFooter
{
	static const uchar DummyByteValue			= 0xCC;
	static const uint32 DatasetTypeSize			= 1 + 1;
	static const uint32 CompressionSettingsSize = 1 + 1 + 1 + 8;

	uchar dummyByte;

	fq::FastqDatasetType datasetType;
	CompressionSettings compSettings;

	enum DatasetTypeFlags
	{
		FLAG_PLUS_REPETITION	= BIT(0),
		FLAG_COLOR_SPACE		= BIT(1)
	};

	enum CompressionFlags
	{
		FLAG_LOSSY_QUALITY		= BIT(0),
		FLAG_CALCULATE_CRC32	= BIT(1)
	};

	std::vector<uint32> blockSizes;

	// TODO: serializer/deserializer
};


class DsrcFileWriter
{
	core::FileStreamWriterExt* fileStream;
	DsrcFileHeader fileHeader;
	DsrcFileFooter fileFooter;

	uint64 currentBlockId;

	fq::StreamsInfo fastqStreamInfo;
	fq::StreamsInfo dsrcStreamInfo;

	void WriteFileHeader();
	void WriteFileFooter();

public:
	DsrcFileWriter();
	~DsrcFileWriter();

	void StartCompress(const std::string& filename_);
	void SetDatasetType(const fq::FastqDatasetType& typeInfo_)
	{
		fileFooter.datasetType = typeInfo_;
	}

	void SetCompressionSettings(const CompressionSettings& settings_)
	{
		fileFooter.compSettings = settings_;
	}

	void WriteNextChunk(const DsrcDataChunk* block_);
	void FinishCompress();

	const fq::StreamsInfo& GetFastqStreamInfo() const
	{
		return fastqStreamInfo;
	}

	const fq::StreamsInfo& GetDsrcStreamInfo() const
	{
		return dsrcStreamInfo;
	}
};


class DsrcFileReader
{
	core::FileStreamReaderExt*	fileStream;
	DsrcFileHeader fileHeader;
	DsrcFileFooter fileFooter;

	uint64 currentBlockId;

	void ReadFileHeader();
	void ReadFileFooter();

public:
	DsrcFileReader();
	~DsrcFileReader();

	void StartDecompress(const std::string& fileName_);

	const fq::FastqDatasetType& GetDatasetType() const
	{
		return fileFooter.datasetType;
	}

	const CompressionSettings& GetCompressionSettings() const
	{
		return fileFooter.compSettings;
	}

	bool ReadNextChunk(DsrcDataChunk* block_);
	void FinishDecompress();

	uint64 BlockCount() const
	{
		return fileHeader.blockCount;
	}
};

} // namespace comp

} // namespace dsrc

#endif
