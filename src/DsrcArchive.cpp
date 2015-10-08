/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include <algorithm>

#include "../include/dsrc/DsrcArchive.h"

#include "RecordsBlockCompressor.h"
#include "FastqParser.h"
#include "DsrcFile.h"
#include "DsrcIo.h"
#include "utils.h"

namespace dsrc
{

using namespace comp;
using namespace core;

namespace ext
{

struct ArchiveImplBase
{
	BlockCompressor* compressor;
	fq::FastqDataChunk* fastqChunk;
	DsrcDataChunk* dsrcChunk;

	CompressionSettings compSettings;

	ArchiveImplBase()
		:	compressor(NULL)
		,	fastqChunk(NULL)
		,	dsrcChunk(NULL)
	{}

	virtual ~ArchiveImplBase()
	{}
};


// DSRC archive writer/reader implementation
//
struct DsrcArchiveWriter::ArchiveWriterImpl : public ArchiveImplBase
{
	DsrcFileWriter* dsrcWriter;
	std::string errorMsg;

	bool isDatasetTypeKnown;
	uint32 fastqQualityOffset;

	ArchiveWriterImpl()
		:	dsrcWriter(NULL)
		,	isDatasetTypeKnown(false)
	{}

	~ArchiveWriterImpl()
	{
		TFree(fastqChunk);
		TFree(compressor);

		TFree(dsrcWriter);
		TFree(dsrcChunk);
	}

	bool StartCompress(const std::string& dsrcFilename_,
					   const DsrcCompressionSettings& compressionSettings_,
					   uint32 qualityOffset_)
	{
		if (IsError())
			ClearError();

		if (!ValidateInputArguments(dsrcFilename_, compressionSettings_, qualityOffset_))
			return false;

		if (dsrcWriter == NULL)
		{
			dsrcWriter = new DsrcFileWriter();
		}

		// TODO: uniform settings
		compSettings = CompressionSettings::ConvertFrom(compressionSettings_);
		dsrcWriter->SetCompressionSettings(compSettings);

		// we will set the dataset type later, explicitely after analysis
		isDatasetTypeKnown = false;
		fastqQualityOffset = qualityOffset_;

		dsrcWriter->SetDatasetType(FastqDatasetType::Default());
		if(compressor == NULL)
		{
			compressor = new BlockCompressor(FastqDatasetType::Default(), compSettings);
			dsrcChunk = new DsrcDataChunk();
		}
		compressor->Reset();

		try
		{
			dsrcWriter->StartCompress(dsrcFilename_.c_str());
		}
		catch (const DsrcException& e)
		{
			AddError(e.what());
		}

		return !IsError();
	}

	void SetFastqDatasetType(const FastqDatasetType& type_)
	{
		ASSERT(compressor != NULL);
		ASSERT(dsrcWriter != NULL);

		compressor->Reconfigure(type_, compSettings);
		dsrcWriter->SetDatasetType(type_);
	}

	void FinishCompress()
	{
		dsrcWriter->FinishCompress();
	}

	bool ValidateInputArguments(const std::string& dsrcFilename_,
								const DsrcCompressionSettings& compressionSettings_,
								uint32 qualityOffset_)
	{
		if (dsrcFilename_.length() == 0)
			AddError("no input DSRC file specified");

		if (compressionSettings_.dnaCompressionLevel > DsrcCompressionSettings::MaxDnaCompressionLevel)
			AddError("invalid DNA compression mode specified [0-3]\n");

		if (compressionSettings_.qualityCompressionLevel > DsrcCompressionSettings::MaxQualityCompressionLevel)
			AddError("invalid Quality compression mode specified [0-2]\n");

		if ( !(compressionSettings_.fastqBufferSizeMb >= DsrcCompressionSettings::MinFastqBufferSizeMB
			   && compressionSettings_.fastqBufferSizeMb <= DsrcCompressionSettings::MaxFastqBufferSizeMB) )
		{
			AddError("invalid fastq buffer size specified [1-1024] \n");
		}

		if (qualityOffset_ != FastqDatasetType::AutoQualityOffsetSelect
				&& !(qualityOffset_ >= 33 && qualityOffset_ <= 64) )
		{
			AddError("invalid Quality offset mode specified [33, 64]");
		}

		return !IsError();
	}

	bool IsError() const
	{
		return errorMsg.length() > 0;
	}

	void AddError(const std::string& err_)
	{
		errorMsg += "Error: " + err_ + '\n';
	}

	void ClearError()
	{
		errorMsg.clear();
	}
};


struct DsrcArchiveReader::ArchiveReaderImpl : public ArchiveImplBase
{
	DsrcFileReader* dsrcReader;
	std::string errorMsg;

	ArchiveReaderImpl()
		:	dsrcReader(NULL)
	{}

	~ArchiveReaderImpl()
	{
		TFree(fastqChunk);
		TFree(compressor);

		TFree(dsrcReader);
		TFree(dsrcChunk);
	}

	bool StartDecompress(const std::string& dsrcFilename_)
	{
		if (IsError())
			ClearError();

		if (!ValidateInputArguments(dsrcFilename_))
			return false;

		if (dsrcReader == NULL)
		{
			dsrcReader = new DsrcFileReader();
		}

		try
		{
			dsrcReader->StartDecompress(dsrcFilename_.c_str());
		}
		catch (const DsrcException& e)
		{
			AddError(e.what());
			return false;
		}

		if (compressor == NULL)
		{
			compSettings = dsrcReader->GetCompressionSettings();
			compressor = new BlockCompressor(dsrcReader->GetDatasetType(),
											 compSettings);

			ASSERT(dsrcChunk == NULL);
			dsrcChunk = new DsrcDataChunk();
		}

		compressor->Reset();

		return true;
	}

	void FinishDecompress()
	{
		dsrcReader->FinishDecompress();
	}

	bool ValidateInputArguments(const std::string& dsrcFilename_)
	{
		if (dsrcFilename_.length() == 0)
			AddError("no input DSRC file specified");

		return !IsError();
	}

	bool IsError() const
	{
		return errorMsg.length() > 0;
	}

	void AddError(const std::string& err_)
	{
		errorMsg += "Error: " + err_ + '\n';
	}

	void ClearError()
	{
		errorMsg.clear();
	}
};


// DSRC archive writer/reader class
//
DsrcArchiveWriter::DsrcArchiveWriter()
	:	archiveImpl(NULL)
{
	archiveImpl = new ArchiveWriterImpl();
}

DsrcArchiveWriter::~DsrcArchiveWriter()
{
	delete archiveImpl;
}

bool DsrcArchiveWriter::IsError() const
{
	return archiveImpl->IsError();
}

const std::string& DsrcArchiveWriter::GetError() const
{
	return archiveImpl->errorMsg;
}

void DsrcArchiveWriter::ClearError()
{
	archiveImpl->ClearError();
}


DsrcArchiveReader::DsrcArchiveReader()
	:	archiveImpl(NULL)
{
	archiveImpl = new ArchiveReaderImpl();
}

DsrcArchiveReader::~DsrcArchiveReader()
{
	delete archiveImpl;
}

bool DsrcArchiveReader::GetCompressionSettings(DsrcCompressionSettings& settings_) const
{
	if (archiveImpl->dsrcReader == NULL)
		return false;

	settings_ = CompressionSettings::ConvertTo(archiveImpl->dsrcReader->GetCompressionSettings());
	return true;
}


bool DsrcArchiveReader::IsError() const
{
	return archiveImpl->IsError();
}

const std::string& DsrcArchiveReader::GetError() const
{
	return archiveImpl->errorMsg;
}

void DsrcArchiveReader::ClearError()
{
	archiveImpl->ClearError();
}


// DSRC archive records writer/reader implementation
//
struct DsrcArchiveRecordsWriter::RecordsWriterImpl
{
	RecordsBlockCompressor* recordsCompressor;
	fq::FastqDataChunk* fastqChunk;

	bool fastqHasPlusRepetition;		// used when parsing

	RecordsWriterImpl()
		:	recordsCompressor(NULL)
		,	fastqChunk(NULL)
		,	fastqHasPlusRepetition(false)
	{}

	~RecordsWriterImpl()
	{
		TFree(fastqChunk);
		TFree(recordsCompressor);
	}

	void CreateCompressorContext(BlockCompressor& compressor_, uint32 fastqBufferSizeMb_)
	{
		ASSERT(fastqChunk == NULL);
		ASSERT(recordsCompressor == NULL);

		if (fastqChunk == NULL)
		{
			fastqChunk = new fq::FastqDataChunk((uint64)fastqBufferSizeMb_ << 20);
			recordsCompressor = new RecordsBlockCompressor(compressor_, *fastqChunk);
		}

		fastqChunk->size = 0;
		recordsCompressor->Reset();
	}

	void FlushChunk(DsrcFileWriter& dsrcWriter_, DsrcDataChunk& dsrcChunk_)
	{
		core::BitMemoryWriter mem(dsrcChunk_.data);
		recordsCompressor->Flush(mem);
		mem.Flush();
		dsrcChunk_.size = mem.Position();
		dsrcWriter_.WriteNextChunk(&dsrcChunk_);
	}
};


struct DsrcArchiveRecordsReader::RecordsReaderImpl
{
	RecordsBlockCompressor* recordsCompressor;
	fq::FastqDataChunk* fastqChunk;

	RecordsReaderImpl()
		:	recordsCompressor(NULL)
		,	fastqChunk(NULL)
	{}

	~RecordsReaderImpl()
	{
		TFree(fastqChunk);
		TFree(recordsCompressor);
	}

	void CreateCompressorContext(BlockCompressor& compressor_, uint32 fastqBufferSizeMb_)
	{
		if (fastqChunk == NULL)
		{
			fastqChunk = new fq::FastqDataChunk((uint64)fastqBufferSizeMb_ << 20);
			recordsCompressor = new RecordsBlockCompressor(compressor_, *fastqChunk);
		}

		fastqChunk->size = 0;
		recordsCompressor->Reset();
	}

	bool FeedChunk(DsrcFileReader& dsrcReader_, DsrcDataChunk& dsrcChunk_)
	{
		if (!dsrcReader_.ReadNextChunk(&dsrcChunk_))
			return false;

		core::BitMemoryReader mem(dsrcChunk_.data.Pointer(), dsrcChunk_.size);
		recordsCompressor->Feed(mem);
		return true;
	}
};


// DSRC archive records writer/reader class
//
DsrcArchiveRecordsWriter::DsrcArchiveRecordsWriter()
{
	writerImpl = new RecordsWriterImpl();
}

DsrcArchiveRecordsWriter::~DsrcArchiveRecordsWriter()
{
	delete writerImpl;
}

bool DsrcArchiveRecordsWriter::StartCompress(const std::string &filename_,
											 const DsrcCompressionSettings &compressionSettings_,
											 uint32 /*threadsNum_*/,
											 uint32 qualityOffset_)
{
	if (!archiveImpl->StartCompress(filename_, compressionSettings_, qualityOffset_))
		return false;

	writerImpl->CreateCompressorContext(*archiveImpl->compressor, compressionSettings_.fastqBufferSizeMb);

	// variables used when analysing FASTQ records
	writerImpl->fastqHasPlusRepetition = false;

	return true;
}

bool DsrcArchiveRecordsWriter::WriteNextRecord(const FastqRecord& rec_)
{
	if (IsError())
		return false;

	const uint64 approxRecordSize = rec_.sequence.length() * 2 + rec_.tag.length() * 2;

	// check the plus repetitions consistency while writing records
	// unfortunately we cannot handle this feature inside FastqParser
	if (rec_.plus.length() > 1)
	{
		if (!writerImpl->fastqHasPlusRepetition)
			writerImpl->fastqHasPlusRepetition = true;
	}
	else if (rec_.plus.length() == 1)
	{
		if (writerImpl->fastqHasPlusRepetition)
		{
			archiveImpl->AddError("inconsistency in \"+\" lines information content across FASTQ data");
			return false;
		}
	}


	// do we have all the data to compresss the next block?
	if (writerImpl->recordsCompressor->RawChunkSize() + approxRecordSize > (uint64)archiveImpl->compSettings.fastqBufferSizeMb << 20)
	{
		// to automatize the dataset type setting, an analysis will be performed after reading the full first block of data
		if (!archiveImpl->isDatasetTypeKnown)
		{
			FastqDatasetType dsType;
			dsType.qualityOffset = archiveImpl->fastqQualityOffset;
			dsType.plusRepetition = writerImpl->fastqHasPlusRepetition;
			if (!writerImpl->recordsCompressor->AnalyzeRecords(dsType.qualityOffset == FastqDatasetType::AutoQualityOffsetSelect,
															   dsType.colorSpace,
															   dsType.qualityOffset))
			{
				archiveImpl->AddError("problem analyzing FASTQ dataset type");
				return false;
			}

			// set the dataset type in DSRC file writer -- was previously in StartCompress()
			archiveImpl->isDatasetTypeKnown = true;
			archiveImpl->SetFastqDatasetType(dsType);
		}

		writerImpl->FlushChunk(*archiveImpl->dsrcWriter, *archiveImpl->dsrcChunk);
	}

	writerImpl->recordsCompressor->WriteNextRecord(rec_);
	return true;
}

void DsrcArchiveRecordsWriter::FinishCompress()
{
	if (writerImpl->recordsCompressor->RawChunkSize() > 0)
	{
		// handle the case of only one block-archive
		// to automatize the dataset type setting, an analysis will be performed after reading the full first block of data
		if (!archiveImpl->isDatasetTypeKnown)
		{
			FastqDatasetType dsType;
			dsType.qualityOffset = archiveImpl->fastqQualityOffset;
			dsType.plusRepetition = writerImpl->fastqHasPlusRepetition;
			if (!writerImpl->recordsCompressor->AnalyzeRecords(dsType.qualityOffset == FastqDatasetType::AutoQualityOffsetSelect,
															   dsType.colorSpace,
															   dsType.qualityOffset))
			{
				archiveImpl->AddError("problem analyzing FASTQ dataset type");

				// we will return here and set an error
				archiveImpl->FinishCompress();
				return;
			}

			// set the dataset type in DSRC file writer -- was previously in StartCompress()
			archiveImpl->isDatasetTypeKnown = true;
			archiveImpl->SetFastqDatasetType(dsType);
		}


		writerImpl->FlushChunk(*archiveImpl->dsrcWriter, *archiveImpl->dsrcChunk);
	}

	archiveImpl->FinishCompress();
}


DsrcArchiveRecordsReader::DsrcArchiveRecordsReader()
{
	readerImpl = new RecordsReaderImpl();
}

DsrcArchiveRecordsReader::~DsrcArchiveRecordsReader()
{
	delete readerImpl;
}

bool DsrcArchiveRecordsReader::StartDecompress(const std::string &filename_,
											   uint32 /*threadsNum_*/)
{
	if (!archiveImpl->StartDecompress(filename_))
		return false;

	readerImpl->CreateCompressorContext(*archiveImpl->compressor, archiveImpl->compSettings.fastqBufferSizeMb);

	return true;
}

void DsrcArchiveRecordsReader::FinishDecompress()
{
	archiveImpl->FinishDecompress();
}

bool DsrcArchiveRecordsReader::ReadNextRecord(FastqRecord& rec_)
{
	if (!readerImpl->recordsCompressor->ReadNextRecord(rec_))
	{
		if (!readerImpl->FeedChunk(*archiveImpl->dsrcReader, *archiveImpl->dsrcChunk))
			return false;
		return readerImpl->recordsCompressor->ReadNextRecord(rec_);
	}
	return true;
}


// DSRC archive block writer/reader implementation
//
struct DsrcArchiveBlocksWriterST::BlockWriterImpl
{
	fq::StreamsInfo rawStreamInfo;
	fq::StreamsInfo compStreamInfo;
	fq::FastqDataChunk* fastqChunk;

	BlockWriterImpl()
		:	fastqChunk(NULL)
	{}

	~BlockWriterImpl()
	{
		TFree(fastqChunk);
	}
};

struct DsrcArchiveBlocksReaderST::BlockReaderImpl
{
	uint64 currentBlockSize;
	uint64 bytesRead;
	fq::FastqDataChunk* fastqChunk;

	BlockReaderImpl()
		:	currentBlockSize(0)
		,	bytesRead(0)
		,	fastqChunk(NULL)
	{}

	~BlockReaderImpl()
	{
		TFree(fastqChunk);
	}
};


// DSRC archive block writer/reader class
//
DsrcArchiveBlocksWriterST::DsrcArchiveBlocksWriterST()
{
	writerImpl = new BlockWriterImpl();
}

DsrcArchiveBlocksWriterST::~DsrcArchiveBlocksWriterST()
{
	delete writerImpl;
}

bool DsrcArchiveBlocksWriterST::StartCompress(const std::string &filename_,
											  const DsrcCompressionSettings &compressionSettings_,
											  uint32 /*threadsNum_*/,
											  uint32 qualityOffset_)
{
	if (!archiveImpl->StartCompress(filename_, compressionSettings_, qualityOffset_))
		return false;

	// WARN: we will create here an additional aux buffer to perform additional
	// copy operations when writing blocks
	//
	// TODO: after validation, we will move to direct on-memory operations when
	// compressing blocks
	//
	if (writerImpl->fastqChunk == NULL)
	{
		writerImpl->fastqChunk = new fq::FastqDataChunk((uint64)compressionSettings_.fastqBufferSizeMb << 20);
	}
	writerImpl->fastqChunk->size = 0;

	return true;
}

void DsrcArchiveBlocksWriterST::FinishCompress()
{
	archiveImpl->FinishCompress();
}

unsigned long DsrcArchiveBlocksWriterST::WriteNextBlock(const char* buffer_, unsigned long bufferSize_)
{
	ASSERT(writerImpl->fastqChunk != NULL);		// TODO: throw exception when unintialized

	ASSERT(buffer_ != NULL);

	if (bufferSize_ == 0)
		return 0;

	// WARN: we will perform here an additional in-memory copy
	//
	// TODO: after validation, we will move to direct in-buffer data compression
	//
	if (writerImpl->fastqChunk->data.Size() < bufferSize_)
	{
		writerImpl->fastqChunk->data.Extend(bufferSize_);		// Shall we add some extension margin?
	}
	std::copy(buffer_, buffer_ + bufferSize_, writerImpl->fastqChunk->data.Pointer());
	writerImpl->fastqChunk->size = bufferSize_;
	//
	// //


	// analyze the FASTQ dataset type and if not done before
	if (!archiveImpl->isDatasetTypeKnown)
	{
		FastqDatasetType datasetType;
		fq::FastqParser parser;
		if (!parser.Analyze(*writerImpl->fastqChunk,
							datasetType,
							archiveImpl->fastqQualityOffset != FastqDatasetType::AutoQualityOffsetSelect))
		{
			archiveImpl->AddError("problem analyzing FASTQ dataset type");
			return 0;
		}

		archiveImpl->SetFastqDatasetType(datasetType);
		archiveImpl->isDatasetTypeKnown = true;
	}


	// compress the FASTQ block
	core::BitMemoryWriter mem(archiveImpl->dsrcChunk->data);
	archiveImpl->compressor->Store(mem,
								   writerImpl->rawStreamInfo,
								   writerImpl->compStreamInfo,
								   *writerImpl->fastqChunk);
	const uint64 compSize = mem.Position();

	archiveImpl->dsrcChunk->size = compSize;
	archiveImpl->dsrcWriter->WriteNextChunk(archiveImpl->dsrcChunk);

	return compSize;
}


DsrcArchiveBlocksReaderST::DsrcArchiveBlocksReaderST()
{
	readerImpl = new BlockReaderImpl();
}

DsrcArchiveBlocksReaderST::~DsrcArchiveBlocksReaderST()
{
	delete readerImpl;
}

bool DsrcArchiveBlocksReaderST::StartDecompress(const std::string &filename_,
												uint32 /*threadsNum*/)
{
	if (!archiveImpl->StartDecompress(filename_))
		return false;

	// WARN: we will create here an additional aux buffer to perform additional
	// copy operations when writing blocks
	//
	// TODO: after validation, we will move to direct on-memory operations when
	// decompressing blocks
	//
	uint64 fastqBufferSizeMb = archiveImpl->dsrcReader->GetCompressionSettings().fastqBufferSizeMb;
	if (readerImpl->fastqChunk == NULL)
	{
		readerImpl->fastqChunk = new fq::FastqDataChunk((uint64)fastqBufferSizeMb << 20);
	}
	readerImpl->fastqChunk->size = 0;

	return true;
}

void DsrcArchiveBlocksReaderST::FinishDecompress()
{
	archiveImpl->FinishDecompress();
}

unsigned long DsrcArchiveBlocksReaderST::ReadNextBlock(char* buffer_, unsigned long bufferSize_)
{
	ASSERT(readerImpl->fastqChunk != NULL);		// TODO: exceptions or '0'

	ASSERT(buffer_ != NULL);

	if (bufferSize_ == 0)
		return 0;

	// WARN: we will perform here an additional in-memory copy
	//
	// TODO: after validation, we will move to direct in-bufer data compression,
	// but firstly we need to expose to the user the size of the FASTQ block
	//
	if (readerImpl->bytesRead == readerImpl->currentBlockSize)
	{
		if (!archiveImpl->dsrcReader->ReadNextChunk(archiveImpl->dsrcChunk))
			return 0;

		// while decompressing the Buffer::Extend() method can be called, when the size of
		// the input buffer is not large enough to store all the data
		core::BitMemoryReader mem(archiveImpl->dsrcChunk->data.Pointer(), archiveImpl->dsrcChunk->size);
		archiveImpl->compressor->Read(mem, *readerImpl->fastqChunk);
		readerImpl->currentBlockSize = readerImpl->fastqChunk->size;
		readerImpl->bytesRead = 0;
	}

	// single-time full copy of buffer
	uint64 bytesCopied = 0;
	if (bufferSize_ > readerImpl->currentBlockSize - readerImpl->bytesRead)
	{
		std::copy(readerImpl->fastqChunk->data.Pointer() + readerImpl->bytesRead,
				  readerImpl->fastqChunk->data.Pointer() + readerImpl->currentBlockSize,
				  buffer_);

		bytesCopied = readerImpl->currentBlockSize - readerImpl->bytesRead;
		readerImpl->bytesRead = readerImpl->currentBlockSize;
	}
	else	// multiple-time copy, chunk by chunk
	{
		std::copy(readerImpl->fastqChunk->data.Pointer() + readerImpl->bytesRead,
				  readerImpl->fastqChunk->data.Pointer() + readerImpl->bytesRead + bufferSize_,
				  buffer_);

		bytesCopied = bufferSize_;
		readerImpl->bytesRead += bufferSize_;
	}

	return bytesCopied;
}

} // namespace ext

} // namespace dsrc
