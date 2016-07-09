/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "DsrcOperator.h"
#include "FastqStream.h"
#include "DsrcFile.h"
#include "DsrcIo.h"
#include "Fastq.h"
#include "FastqIo.h"
#include "DsrcWorker.h"
#include "FastqParser.h"
#include "BlockCompressor.h"
#include "utils.h"
#include "ErrorHandler.h"

#include <stdexcept>
#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
namespace th = boost;
#else
#include <thread>
namespace th = std;
#endif

#include <sstream>
#include <iomanip>

namespace dsrc
{

namespace comp
{

using namespace core;
using namespace fq;

const uint32 IDsrcOperator::AvailableHardwareThreadsNum = th::thread::hardware_concurrency();


bool DsrcCompressorST::Process(const std::string& fastqFilename_,
							   const std::string& dsrcFilename_,
							   const DsrcCompressionSettings& compSettings_,
							   bool useFastqStdIo_,
							   uint32 qualityOffset_)
{
	ASSERT(!IsError());

	IFastqStreamReader* reader = NULL;
	DsrcFileWriter* writer = NULL;

	// make reusable
	//
	DsrcDataChunk* dsrcChunk = NULL;
	FastqDataChunk* fastqChunk = NULL;
	//
	//

	FastqDatasetType datasetType;
	CompressionSettings settings = CompressionSettings::ConvertFrom(compSettings_);

	try
	{
		if (useFastqStdIo_)
			reader = new FastqStdIoReader();
		else
			reader = new FastqFileReader(fastqFilename_);

		// join into constructor for RAII style
		writer = new DsrcFileWriter();
		writer->StartCompress(dsrcFilename_);

		dsrcChunk = new DsrcDataChunk(DsrcDataChunk::DefaultBufferSize);
		fastqChunk = new FastqDataChunk(compSettings_.fastqBufferSizeMb << 20);

		// analyze the header
		//
		const bool findQOffset = qualityOffset_ == FastqDatasetType::AutoQualityOffsetSelect;
		if (!findQOffset)
		{
			datasetType.qualityOffset = qualityOffset_;
		}

		FastqParser parser;
		if (!reader->ReadNextChunk(fastqChunk) || parser.Analyze(*fastqChunk, datasetType, findQOffset) == 0)
		{
			AddError("Error analyzing FASTQ dataset");
		}
		else
		{
			writer->SetDatasetType(datasetType);
			writer->SetCompressionSettings(settings);
		}
	}
	catch (const DsrcException& e_)
	{
		AddError(e_.what());
	}

	if (!IsError())
	{
		BitMemoryWriter bitMemory(dsrcChunk->data);
		BlockCompressor superblock(datasetType, settings);

		do
		{
			superblock.Store(bitMemory, dsrcChunk->rawStreamsInfo, dsrcChunk->compStreamsInfo, *fastqChunk);

			bitMemory.Flush();
			dsrcChunk->size = bitMemory.Position();

			writer->WriteNextChunk(dsrcChunk);

			if (compSettings_.calculateCrc32)
			{
				BitMemoryReader reader(dsrcChunk->data.Pointer(), dsrcChunk->data.Size());
				std::fill(fastqChunk->data.Pointer(), fastqChunk->data.Pointer() + fastqChunk->data.Size(), 0xCC);

				if (!superblock.VerifyChecksum(reader, *fastqChunk))
				{
					AddError("CRC32 checksums mismatch.");
					break;
				}
			}

			fastqChunk->Reset();
			dsrcChunk->Reset();
			bitMemory.Reset();
		}
		while (reader->ReadNextChunk(fastqChunk));

		reader->Close();
		writer->FinishCompress();


		// set log
		//
		fq::StreamsInfo rawSize = writer->GetFastqStreamInfo();
		fq::StreamsInfo compSize = writer->GetDsrcStreamInfo();

		std::ostringstream ss;
		ss << "Compressed streams sizes (in bytes)\n";
		ss << "TAG: " << std::setw(16) << compSize.sizes[fq::StreamsInfo::MetaStream] + compSize.sizes[fq::StreamsInfo::TagStream]
					  << " / " << std::setw(16) << rawSize.sizes[fq::StreamsInfo::TagStream] << '\n';
		ss << "DNA: " << std::setw(16) << compSize.sizes[fq::StreamsInfo::DnaStream]
					  << " / " << std::setw(16) << rawSize.sizes[fq::StreamsInfo::DnaStream] << '\n';
		ss << "QUA: " << std::setw(16) << compSize.sizes[fq::StreamsInfo::QualityStream]
					  << " / " << std::setw(16) << rawSize.sizes[fq::StreamsInfo::QualityStream] << '\n';
		AddLog(ss.str());
	}

	// make reusable
	//
	TFree(fastqChunk);
	TFree(dsrcChunk);
	//
	//

	TFree(writer);
	TFree(reader);

	return !IsError();
}

bool DsrcDecompressorST::Process(const std::string& fastqFilename_,
								 const std::string& dsrcFilename_,
								 bool useFastqStdIo_)
{
	ASSERT(!IsError());

	DsrcFileReader* reader = NULL;
	IFastqStreamWriter* writer = NULL;

	// make reusable
	//
	DsrcDataChunk* dsrcChunk = NULL;
	FastqDataChunk* fastqChunk = NULL;
	//
	//

	try
	{
		reader = new DsrcFileReader();
		reader->StartDecompress(dsrcFilename_);

		if (useFastqStdIo_)
			writer = new FastqStdIoWriter();
		else
			writer = new FastqFileWriter(fastqFilename_);

		dsrcChunk = new DsrcDataChunk(DsrcDataChunk::DefaultBufferSize);
		fastqChunk = new FastqDataChunk(FastqDataChunk::DefaultBufferSize);
	}
	catch (const DsrcException& e_)
	{
		AddError(e_.what());
	}

	if (!IsError())
	{
		BlockCompressor superblock(reader->GetDatasetType(), reader->GetCompressionSettings());

		while (reader->ReadNextChunk(dsrcChunk))
		{
			BitMemoryReader bitMemory(dsrcChunk->data.Pointer(), dsrcChunk->size);

			superblock.Read(bitMemory, *fastqChunk);

			writer->WriteNextChunk(fastqChunk);

			fastqChunk->Reset();
			dsrcChunk->Reset();
		}

		reader->FinishDecompress();
		writer->Close();
	}

	// make reusable
	//
	TFree(fastqChunk);
	TFree(dsrcChunk);
	//
	//

	TFree(writer);
	TFree(reader);

	return !IsError();
}

bool DsrcCompressorMT::Process(const std::string& fastqFilename_,
							   const std::string& dsrcFilename_,
							   const DsrcCompressionSettings& compSettings_,
							   uint32 threadNum_,
							   bool useFastqStdIo_,
							   uint32 qualityOffset_)
{
	IFastqStreamReader* fileReader = NULL;
	DsrcFileWriter* fileWriter = NULL;

	// make reusable
	//
	FastqDataPool* fastqPool = NULL;
	FastqDataQueue* fastqQueue = NULL;
	DsrcDataPool* dsrcPool = NULL;
	DsrcDataQueue* dsrcQueue = NULL;
	ErrorHandler* errorHandler = NULL;
	//
	//

	FastqReader* dataReader = NULL;
	DsrcWriter* dataWriter = NULL;

	FastqDatasetType datasetType;
	CompressionSettings compSettings = CompressionSettings::ConvertFrom(compSettings_);

	try
	{
		if (useFastqStdIo_)
			fileReader = new FastqStdIoReader();
		else
			fileReader = new FastqFileReader(fastqFilename_);

		fileWriter = new DsrcFileWriter();
		fileWriter->StartCompress(dsrcFilename_);

		const uint32 partNum = (compSettings_.fastqBufferSizeMb < 128) ? threadNum_ * 4 : threadNum_ * 2;
		fastqPool = new FastqDataPool(partNum, compSettings_.fastqBufferSizeMb << 20);		// maxPart, bufferPartSize
		fastqQueue = new FastqDataQueue(partNum, 1);										// maxPart, threadCount

		dsrcPool = new DsrcDataPool(partNum, compSettings_.fastqBufferSizeMb << 20);
		dsrcQueue = new DsrcDataQueue(partNum, threadNum_);

		if (compSettings_.calculateCrc32)
			errorHandler = new MultithreadedErrorHandler();
		else
			errorHandler = new ErrorHandler();

		dataReader = new FastqReader(*fileReader, *fastqQueue, *fastqPool, *errorHandler);
		dataWriter = new DsrcWriter(*fileWriter, *dsrcQueue, *dsrcPool, *errorHandler);

		// analyze file
		//
		const bool findQOffset = qualityOffset_ == FastqDatasetType::AutoQualityOffsetSelect;
		if (!findQOffset)
			datasetType.qualityOffset = qualityOffset_;

		if (!dataReader->AnalyzeFirstChunk(datasetType, findQOffset))
		{
			AddError("Error analyzing FASTQ dataset");
		}
		else
		{
			fileWriter->SetDatasetType(datasetType);
			fileWriter->SetCompressionSettings(compSettings);
		}
	}
	catch (const DsrcException& e_)
	{
		AddError(e_.what());
	}

	if (!IsError())
	{
		// launch threads
		//
		th::thread readerThread(th::ref(*dataReader));

		std::vector<DsrcCompressor*> operators;
		operators.resize(threadNum_);

#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;		// why C++11 does not have thread_group? ://

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new DsrcCompressor(*fastqQueue, *fastqPool, *dsrcQueue, *dsrcPool, *errorHandler, datasetType, compSettings);
			opThreadGroup.create_thread(th::ref(*operators[i]));
		}

		(*dataWriter)();			// main thread works as writer

		readerThread.join();
		opThreadGroup.join_all();

#else
		std::vector<th::thread> opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new DsrcCompressor(*fastqQueue, *fastqPool, *dsrcQueue, *dsrcPool, *errorHandler, datasetType, compSettings);
			opThreadGroup.push_back(th::thread(th::ref(*operators[i])));
		}

		(*dataWriter)();

		readerThread.join();

		// find difference: 'for (std::vector<th::thread>::iterator i = opThreadGroup.begin(); i != opThreadGroup.end(); ++i)'
		// --> 'for (auto& t : opThreadGroup)'
		for (th::thread& t : opThreadGroup)
		{
			t.join();
		}

#endif

		// check for errors
		//
		if (errorHandler->IsError())
			AddError(errorHandler->GetError());


		// free resources, cleanup
		//
		fastqQueue->Reset();
		dsrcQueue->Reset();

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			delete operators[i];
		}

		fileReader->Close();
		fileWriter->FinishCompress();


		// set log
		//
		fq::StreamsInfo rawSize = fileWriter->GetFastqStreamInfo();
		fq::StreamsInfo compSize = fileWriter->GetDsrcStreamInfo();

		std::ostringstream ss;
		ss << "Compressed streams sizes (in bytes)\n";
		ss << "TAG: " << std::setw(16) << compSize.sizes[fq::StreamsInfo::MetaStream] + compSize.sizes[fq::StreamsInfo::TagStream]
					  << " / " << std::setw(16) << rawSize.sizes[fq::StreamsInfo::TagStream] << '\n';
		ss << "DNA: " << std::setw(16) << compSize.sizes[fq::StreamsInfo::DnaStream]
					  << " / " << std::setw(16) << rawSize.sizes[fq::StreamsInfo::DnaStream] << '\n';
		ss << "QUA: " << std::setw(16) << compSize.sizes[fq::StreamsInfo::QualityStream]
					  << " / " << std::setw(16) << rawSize.sizes[fq::StreamsInfo::QualityStream] << '\n';
		AddLog(ss.str());
	}

	TFree(dataWriter);
	TFree(dataReader);
	TFree(errorHandler);

	// make reusable
	//
	TFree(dsrcQueue);
	TFree(dsrcPool);
	TFree(fastqQueue);
	TFree(fastqPool);
	//
	//

	TFree(fileWriter);
	TFree(fileReader);

	return !IsError();
}

bool DsrcDecompressorMT::Process(const std::string& fastqFilename_,
								 const std::string& dsrcFilename_,
								 uint32 threadNum_,
								 bool useFastqStdIo_)
{
	ASSERT(!IsError());

	DsrcFileReader* fileReader = NULL;
	IFastqStreamWriter* fileWriter = NULL;

	// make reusable
	//
	FastqDataPool* fastqPool = NULL;
	FastqDataQueue* fastqQueue = NULL;
	DsrcDataPool* dsrcPool = NULL;
	DsrcDataQueue* dsrcQueue = NULL;
	ErrorHandler* errorHandler = NULL;
	//
	//

	DsrcReader* dataReader = NULL;
	FastqWriter* dataWriter = NULL;

	try
	{
		fileReader = new DsrcFileReader();
		fileReader->StartDecompress(dsrcFilename_);		// here get FASTQ buffer size!!!!

		if (useFastqStdIo_)
			fileWriter= new FastqStdIoWriter();
		else
			fileWriter = new FastqFileWriter(fastqFilename_);

		const uint32 fastqBufferSizeMB = fileReader->GetCompressionSettings().fastqBufferSizeMb;
		const uint32 partNum = (fastqBufferSizeMB < 128) ? threadNum_ * 4 : threadNum_ * 2;

		dsrcPool = new DsrcDataPool(partNum, fastqBufferSizeMB << 20);
		dsrcQueue = new DsrcDataQueue(partNum, 1);

		fastqPool = new FastqDataPool(partNum, fastqBufferSizeMB << 20);		// maxPart, bufferPartSize
		fastqQueue = new FastqDataQueue(partNum, threadNum_);					// maxPart, threadCount

		errorHandler = new ErrorHandler();
		dataReader = new DsrcReader(*fileReader, *dsrcQueue, *dsrcPool, *errorHandler);
		dataWriter = new FastqWriter(*fileWriter, *fastqQueue, *fastqPool, *errorHandler);
	}
	catch (const DsrcException& e_)
	{
		AddError(e_.what());
	}

	if (!IsError())
	{
		// launch threads
		//
		th::thread readerThread(th::ref(*dataReader));

		std::vector<DsrcDecompressor*> operators;
		operators.resize(threadNum_);


#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new DsrcDecompressor(*fastqQueue, *fastqPool, *dsrcQueue, *dsrcPool, *errorHandler,
												fileReader->GetDatasetType(), fileReader->GetCompressionSettings());
			opThreadGroup.create_thread(th::ref(*operators[i]));
		}

		(*dataWriter)();	// main thread works as writer

		readerThread.join();

		opThreadGroup.join_all();
#else
		std::vector<th::thread> opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new DsrcDecompressor(*fastqQueue, *fastqPool, *dsrcQueue, *dsrcPool, *errorHandler,
												fileReader->GetDatasetType(), fileReader->GetCompressionSettings());
			opThreadGroup.push_back(th::thread(th::ref(*operators[i])));
		}

		(*dataWriter)();

		readerThread.join();

		for (th::thread& t : opThreadGroup)
		{
			t.join();
		}
#endif

		// free resources, cleanup
		//
		fastqQueue->Reset();
		dsrcQueue->Reset();

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			delete operators[i];
		}

		fileReader->FinishDecompress();
		fileWriter->Close();
	}

	TFree(dataWriter);
	TFree(dataReader);
	TFree(errorHandler);

	// make reusable
	//
	TFree(dsrcQueue);
	TFree(dsrcPool);
	TFree(fastqQueue);
	TFree(fastqPool);
	//
	//

	TFree(fileWriter);
	TFree(fileReader);

	return !IsError();
}

} // namespace comp

} // namespace dsrc
