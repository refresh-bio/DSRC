/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "../include/dsrc/Globals.h"

#include <iostream>

#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
namespace th = boost;
#else
#include <thread>
namespace th = std;
#endif

#include "../include/dsrc/DsrcModule.h"

#include "DsrcOperator.h"


namespace dsrc
{

using namespace comp;

namespace ext
{

const uint32 DsrcModule::AvailableHardwareThreadsNum = th::thread::hardware_concurrency();



// TODO: externalize the version and build date string in makefile
std::string DsrcModule::Version()
{
	return std::string("2.1.0 @ 05.10.2015");
}


DsrcModule::DsrcModule()
{}

DsrcModule::~DsrcModule()
{}

bool DsrcModule::Compress(const std::string &inFastqFilename_,
						  const std::string &outDsrcFilename_,
						  const DsrcCompressionSettings& compSettings_,
						  uint32 threadsNum_,
						  bool useFastqStdIo_,
						  uint32 qualityOffset_)
{
	if (IsError())
		ClearError();

	// check input params
	//
	if (outDsrcFilename_.length() == 0)
		AddError("no output DSRC file specified");

	if (inFastqFilename_.length() == 0 && !useFastqStdIo_)
		AddError("no input FASTQ file specified");

	if (compSettings_.dnaCompressionLevel > DsrcCompressionSettings::MaxDnaCompressionLevel)
		AddError("invalid DNA compression mode specified [0-3]\n");

	if (compSettings_.qualityCompressionLevel > DsrcCompressionSettings::MaxQualityCompressionLevel)
		AddError("invalid Quality compression mode specified [0-2]\n");

	if ( !(compSettings_.fastqBufferSizeMb >= DsrcCompressionSettings::MinFastqBufferSizeMB
		   && compSettings_.fastqBufferSizeMb <= DsrcCompressionSettings::MaxFastqBufferSizeMB) )
	{
		AddError("invalid fastq buffer size specified [1-1024] \n");
	}

	if (qualityOffset_ != FastqDatasetType::AutoQualityOffsetSelect
			&& !(qualityOffset_ >= 33 && qualityOffset_ <= 64) )
	{
		AddError("invalid Quality offset mode specified [33, 64]");
	}

	if (IsError())
		return false;


	// compress
	//
	if (threadsNum_ > 0)
	{
		DsrcCompressorMT dsrc;
		if (!dsrc.Process(inFastqFilename_, outDsrcFilename_, compSettings_,
						  threadsNum_, useFastqStdIo_, qualityOffset_))
		{
			SetError(dsrc.GetError());
			return false;
		}
		if (dsrc.GetLog().length() > 0)
			AddLog(dsrc.GetLog());
	}
	else
	{
		DsrcCompressorST dsrc;
		if (!dsrc.Process(inFastqFilename_, outDsrcFilename_, compSettings_,
						  useFastqStdIo_, qualityOffset_))
		{
			SetError(dsrc.GetError());
			return false;
		}
		if (dsrc.GetLog().length() > 0)
			AddLog(dsrc.GetLog());
	}
	return true;
}

bool DsrcModule::Decompress(const std::string &inDsrcFilename_,
							const std::string &outFastqFilename_,
							uint32 threadsNum_,
							bool useFastqStdIo_)
{
	if (IsError())
		ClearError();

	// check input params
	//
	if (inDsrcFilename_.length() == 0)
		AddError("no input DSRC file specified");

	if (outFastqFilename_.length() == 0 && !useFastqStdIo_)
		AddError("no input FASTQ file specified");

	if (IsError())
		return false;


	// decompress
	//
	if (threadsNum_ > 0)
	{
		DsrcDecompressorMT dsrc;
		if (!dsrc.Process(outFastqFilename_,
						  inDsrcFilename_,
						  threadsNum_,
						  useFastqStdIo_))
		{
			SetError(dsrc.GetError());
			return false;
		}
		if (dsrc.GetLog().length() > 0)
			AddLog(dsrc.GetLog());
	}
	else
	{
		DsrcDecompressorST dsrc;
		if (!dsrc.Process(outFastqFilename_,
						  inDsrcFilename_,
						  useFastqStdIo_))
		{
			SetError(dsrc.GetError());
			return false;
		}
		if (dsrc.GetLog().length() > 0)
			AddLog(dsrc.GetLog());
	}
	return true;
}


// error handling
//
bool DsrcModule::IsError() const
{
	return errorMsg.length() > 0;
}

const std::string& DsrcModule::GetError() const
{
	return errorMsg;
}

void DsrcModule::AddError(const std::string& err_)
{
	ASSERT(err_.length() > 0);
	errorMsg += "Error: " + err_ + "\n";
}

void DsrcModule::SetError(const std::string& err_)
{
	ASSERT(err_.length() > 0);
	errorMsg = err_;
}

void DsrcModule::ClearError()
{
	errorMsg.clear();
}

const std::string& DsrcModule::GetLog() const
{
	return logMsg;
}

void DsrcModule::AddLog(const std::string& log_)
{
	ASSERT(log_.length() > 0);
	logMsg += log_ + "\n";
}

void DsrcModule::ClearLog()
{
	logMsg.clear();
}

} // namespace ext

} // namespace dsrc
