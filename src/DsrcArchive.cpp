/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "../include/dsrc/DsrcArchive.h"

#include "BlockCompressorExt.h"
#include "DsrcFile.h"
#include "DsrcIo.h"

namespace dsrc
{

using namespace comp;

namespace wrap
{

struct ArchiveSettings
{
	CompressionSettings compSettings;
	fq::FastqDatasetType fastqSettings;
	uint64 fastqBufferSize;

	ArchiveSettings()
		:	fastqBufferSize(0)
	{}

	void FromInputParams(const Configurable& params_)
	{
		ASSERT(compSettings.dnaOrder <= 3);
		ASSERT(compSettings.qualityOrder == 0 || (compSettings.lossy && compSettings.qualityOrder <= 2));
		ASSERT(fastqSettings.qualityOffset == 33 || fastqSettings.qualityOffset == 64);
		ASSERT(params_.GetFastqBufferSizeMB() > 0 && params_.GetFastqBufferSizeMB() < 1024);

		compSettings.dnaOrder = params_.GetDnaCompressionLevel() * 3;
		compSettings.qualityOrder = params_.GetQualityCompressionLevel() * 3;
		compSettings.lossy = params_.IsLossyCompression();

		fastqSettings.qualityOffset = params_.GetQualityOffset();
		fastqSettings.plusRepetition = params_.IsPlusRepetition();
		fastqSettings.colorSpace = params_.IsColorSpace();

		fastqBufferSize = params_.GetFastqBufferSizeMB() << 20UL;
	}

	void ToInputParams(Configurable& params_)
	{
		ASSERT(compSettings.dnaOrder <= 9);
		ASSERT(compSettings.qualityOrder == 0 || (compSettings.lossy && compSettings.qualityOrder <= 6));
		ASSERT(fastqSettings.qualityOffset == 33 || fastqSettings.qualityOffset == 64);
		ASSERT((fastqBufferSize >> 20UL) > 0 && (fastqBufferSize >> 20UL) < 1024UL);

		params_.SetDnaCompressionLevel(compSettings.dnaOrder / 3);
		params_.SetQualityCompressionLevel(compSettings.qualityOrder / 3);
		params_.SetLossyCompression(compSettings.lossy);

		params_.SetQualityOffset(fastqSettings.qualityOffset);
		params_.SetPlusRepetition(fastqSettings.plusRepetition);
		params_.SetColorSpace(fastqSettings.colorSpace);

		params_.SetFastqBufferSizeMB(fastqBufferSize >> 20UL);
	}
};


struct DsrcArchive::ArchiveImpl
{
	BlockCompressorExt* compressor;
	DsrcFileReader* dsrcReader;
	DsrcFileWriter* dsrcWriter;
	DsrcDataChunk* dsrcChunk;

	ArchiveSettings settings;

	ArchiveImpl()
		:	compressor(NULL)
		,	dsrcReader(NULL)
		,	dsrcWriter(NULL)
		,	dsrcChunk(NULL)
	{}
};

DsrcArchive::DsrcArchive()
	:	state(StateNone)
	,	impl(NULL)
{
	impl = new ArchiveImpl();

	// set default parameters
	impl->settings.fastqBufferSize = InputParameters::DefaultFastqBufferSizeMB << 20;

	impl->settings.compSettings = CompressionSettings::Default();
	impl->settings.fastqSettings = fq::FastqDatasetType::Default();
}

DsrcArchive::~DsrcArchive()
{
	if (impl->dsrcReader != NULL)
		delete impl->dsrcReader;
	if (impl->dsrcWriter != NULL)
		delete impl->dsrcWriter;
	if (impl->compressor != NULL)
		delete impl->compressor;
	if (impl->dsrcChunk != NULL)
		delete impl->dsrcChunk;

	delete impl;
}

void DsrcArchive::StartCompress(const std::string& filename_)
{
	ASSERT(state == StateNone);

	impl->settings.FromInputParams(*this);

	if (impl->dsrcWriter == NULL)
	{	
		impl->dsrcWriter = new DsrcFileWriter();
	}
	impl->dsrcWriter->SetCompressionSettings(impl->settings.compSettings);
	impl->dsrcWriter->SetDatasetType(impl->settings.fastqSettings);

	if(impl->compressor == NULL)
	{
		impl->compressor = new BlockCompressorExt(impl->settings.fastqSettings, impl->settings.compSettings);

		ASSERT(impl->dsrcChunk == NULL);
		impl->dsrcChunk = new DsrcDataChunk();
	}

	impl->compressor->Reset();

	impl->dsrcWriter->StartCompress(filename_.c_str());

	state = StateCompression;
}

void DsrcArchive::WriteNextRecord(const FastqRecord& rec_)
{
	ASSERT(state == StateCompression);

	impl->compressor->WriteNextRecord(rec_);

	if (impl->compressor->ChunkSize() > impl->settings.fastqBufferSize)
	{
		FlushChunk();
	}
}

void DsrcArchive::FinishCompress()
{
	ASSERT(state == StateCompression);

	if (impl->compressor->ChunkSize() > 0)
	{
		FlushChunk();
	}

	impl->dsrcWriter->FinishCompress();

	state = StateNone;
}

void DsrcArchive::StartDecompress(const std::string& filename_)
{
	ASSERT(state == StateNone);

	if (impl->dsrcReader == NULL)
	{
		impl->dsrcReader = new DsrcFileReader();
	}
	impl->dsrcReader->StartDecompress(filename_.c_str());

	impl->settings.fastqSettings = impl->dsrcReader->GetDatasetType();
	impl->settings.compSettings = impl->dsrcReader->GetCompressionSettings();
	impl->settings.ToInputParams(*this);

	if (impl->compressor == NULL)
	{
		impl->compressor = new BlockCompressorExt(impl->settings.fastqSettings, impl->settings.compSettings);

		ASSERT(impl->dsrcChunk == NULL);
		impl->dsrcChunk = new DsrcDataChunk();
	}

	impl->compressor->Reset();

	state = StateDecompression;
}

bool DsrcArchive::ReadNextRecord(FastqRecord& rec_)
{
	ASSERT(state == StateDecompression);

	if (!impl->compressor->ReadNextRecord(rec_))
	{
		if (!FeedChunk())
			return false;
		return impl->compressor->ReadNextRecord(rec_);
	}
	return true;
}

void DsrcArchive::FinishDecompress()
{
	ASSERT(state == StateDecompression);

	impl->dsrcReader->FinishDecompress();
}

void DsrcArchive::FlushChunk()
{
	core::BitMemoryWriter mem(impl->dsrcChunk->data);
	impl->compressor->Flush(mem);
	mem.Flush();
	impl->dsrcChunk->size = mem.Position();
	impl->dsrcWriter->WriteNextChunk(impl->dsrcChunk);
}

bool DsrcArchive::FeedChunk()
{
	if (!impl->dsrcReader->ReadNextChunk(impl->dsrcChunk))
		return false;

	core::BitMemoryReader mem(impl->dsrcChunk->data.Pointer(), impl->dsrcChunk->size);
	impl->compressor->Feed(mem);
	return true;
}

} // namespace wrap

} // namespace dsrc
