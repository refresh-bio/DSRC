/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCOPERATOR
#define H_DSRCOPERATOR

#include "../include/dsrc/Globals.h"

#include <string>

#include "FastqStream.h"
#include "DsrcFile.h"

namespace dsrc
{

namespace comp
{

class IDsrcOperator
{
public:
	static const uint32 AvailableHardwareThreadsNum;

	virtual ~IDsrcOperator() {}

	bool IsError() const
	{
		return errorMessage.length() > 0;
	}

	const std::string& GetError() const
	{
		return errorMessage;
	}

	void ClearError()
	{
		errorMessage.clear();
	}

	const std::string& GetLog() const
	{
		return logMessage;
	}

	void ClearLog()
	{
		logMessage.clear();
	}

protected:
	std::string errorMessage;
	std::string logMessage;

	void AddError(const std::string& err_)
	{
		errorMessage += "Error: " + err_ + '\n';
	}

	void AddLog(const std::string& log_)
	{
		logMessage += log_ + '\n';
	}
};

class DsrcCompressorST : public IDsrcOperator
{
public:
	bool Process(const std::string& fastqFilename_,
				 const std::string& dsrcFilename_,
				 const DsrcCompressionSettings& compSettings_,
				 bool useFastqStdIo_,
				 uint32 qualityOffset_);
};

class DsrcDecompressorST : public IDsrcOperator
{
public:
	bool Process(const std::string& fastqFilename_,
				 const std::string& dsrcFilename_,
				 bool useFastqStdIo_);
};

class DsrcCompressorMT : public IDsrcOperator
{
public:
	bool Process(const std::string& fastqFilename_,
				 const std::string& dsrcFilename_,
				 const DsrcCompressionSettings& compSettings_,
				 uint32 threadNum_,
				 bool useFastqStdIo_,
				 uint32 qualityOffset_);
};

class DsrcDecompressorMT : public IDsrcOperator
{
public:
	bool Process(const std::string& fastqFilename_,
				 const std::string& dsrcFilename_,
				 uint32 threadNum_,
				 bool useFastqStdIo_);
};

} // namespace comp

} // namespace dsrc

#endif
