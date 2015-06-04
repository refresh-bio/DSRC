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


bool DsrcCompressorST::Process(const InputParameters &args_)
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
	CompressionSettings settings = GetCompressionSettings(args_);

	try
	{
		if (args_.useFastqStdIo)
			reader = new FastqStdIoReader();
		else
			reader = new FastqFileReader(args_.inputFilename);

		// join into constructor for RAII style
		writer = new DsrcFileWriter();
		writer->StartCompress(args_.outputFilename);

		dsrcChunk = new DsrcDataChunk(DsrcDataChunk::DefaultBufferSize);
		fastqChunk = new FastqDataChunk(args_.fastqBufferSizeMB << 20);

		// analyze the header
		//
		const bool findQOffset = args_.qualityOffset == fq::FastqDatasetType::AutoQualityOffset;
		if (!findQOffset)
		{
			datasetType.qualityOffset = args_.qualityOffset;
		}

		FastqParser parser;
		if (!reader->ReadNextChunk(fastqChunk) || parser.Analyze(*fastqChunk, datasetType, findQOffset) == 0)
		{
			throw DsrcException("Error analyzing FASTQ dataset");
		}

		writer->SetDatasetType(datasetType);
		writer->SetCompressionSettings(settings);
	}
	catch (const std::runtime_error& e_)
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

			if (args_.calculateCrc32)
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

bool DsrcDecompressorST::Process(const InputParameters& args_)
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
		reader->StartDecompress(args_.inputFilename);

		if (args_.useFastqStdIo)
			writer = new FastqStdIoWriter();
		else
			writer = new FastqFileWriter(args_.outputFilename);

		dsrcChunk = new DsrcDataChunk(DsrcDataChunk::DefaultBufferSize);
		fastqChunk = new FastqDataChunk(FastqDataChunk::DefaultBufferSize);
	}
	catch (const std::runtime_error& e_)
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

bool DsrcCompressorMT::Process(const InputParameters &args_)
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
	CompressionSettings compSettings = GetCompressionSettings(args_);

	try
	{
		if (args_.useFastqStdIo)
			fileReader = new FastqStdIoReader();
		else
			fileReader = new FastqFileReader(args_.inputFilename);

		fileWriter = new DsrcFileWriter();
		fileWriter->StartCompress(args_.outputFilename);

		const uint32 partNum = (args_.fastqBufferSizeMB < 128) ? args_.threadNum * 4 : args_.threadNum * 2;
		fastqPool = new FastqDataPool(partNum, args_.fastqBufferSizeMB << 20);		// maxPart, bufferPartSize
		fastqQueue = new FastqDataQueue(partNum, 1);								// maxPart, threadCount

		dsrcPool = new DsrcDataPool(partNum, args_.fastqBufferSizeMB << 20);
		dsrcQueue = new DsrcDataQueue(partNum, args_.threadNum);

		if (args_.calculateCrc32)
			errorHandler = new MultithreadedErrorHandler();
		else
			errorHandler = new ErrorHandler();

		dataReader = new FastqReader(*fileReader, *fastqQueue, *fastqPool, *errorHandler);
		dataWriter = new DsrcWriter(*fileWriter, *dsrcQueue, *dsrcPool, *errorHandler);

		// analyze file
		//
		const bool findQOffset = args_.qualityOffset == fq::FastqDatasetType::AutoQualityOffset;
		if (!findQOffset)
			datasetType.qualityOffset = args_.qualityOffset;

		if (!dataReader->AnalyzeFirstChunk(datasetType, findQOffset))
		{
			throw DsrcException("Error analyzing FASTQ dataset");
		}

		fileWriter->SetDatasetType(datasetType);
		fileWriter->SetCompressionSettings(compSettings);
	}
	catch (const std::runtime_error& e_)
	{
		AddError(e_.what());
	}

	if (!IsError())
	{
		const uint32 threadsNum = args_.threadNum;

		// launch threads
		//
		th::thread readerThread(th::ref(*dataReader));

		std::vector<DsrcCompressor*> operators;
		operators.resize(threadsNum);

#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;		// why C++11 does not have thread_group? ://

		for (uint32 i = 0; i < threadsNum; ++i)
		{
			operators[i] = new DsrcCompressor(*fastqQueue, *fastqPool, *dsrcQueue, *dsrcPool, *errorHandler, datasetType, compSettings);
			opThreadGroup.create_thread(th::ref(*operators[i]));
		}

		(*dataWriter)();			// main thread works as writer

		readerThread.join();
		opThreadGroup.join_all();

#else
		std::vector<th::thread> opThreadGroup;

		for (uint32 i = 0; i < threadsNum; ++i)
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

		for (uint32 i = 0; i < threadsNum; ++i)
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

bool DsrcDecompressorMT::Process(const InputParameters &args_)
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
		fileReader->StartDecompress(args_.inputFilename);

		if (args_.useFastqStdIo)
			fileWriter= new FastqStdIoWriter();
		else
			fileWriter = new FastqFileWriter(args_.outputFilename);

		const uint32 partNum = (args_.fastqBufferSizeMB < 128) ? args_.threadNum * 4 : args_.threadNum * 2;
		dsrcPool = new DsrcDataPool(partNum, args_.fastqBufferSizeMB << 20);
		dsrcQueue = new DsrcDataQueue(partNum, 1);

		fastqPool = new FastqDataPool(partNum, DsrcDataPool::DefaultBufferPartSize);		// maxPart, bufferPartSize
		fastqQueue = new FastqDataQueue(partNum, args_.threadNum);						// maxPart, threadCount

		errorHandler = new ErrorHandler();
		dataReader = new DsrcReader(*fileReader, *dsrcQueue, *dsrcPool, *errorHandler);
		dataWriter = new FastqWriter(*fileWriter, *fastqQueue, *fastqPool, *errorHandler);
	}
	catch(const std::runtime_error& e_)
	{
		AddError(e_.what());
	}

	if (!IsError())
	{
		const uint32 threadsNum = args_.threadNum;

		// launch threads
		//
		th::thread readerThread(th::ref(*dataReader));

		std::vector<DsrcDecompressor*> operators;
		operators.resize(threadsNum);


#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;

		for (uint32 i = 0; i < threadsNum; ++i)
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

		for (uint32 i = 0; i < threadsNum; ++i)
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

		for (uint32 i = 0; i < threadsNum; ++i)
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
