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

	virtual bool Process(const InputParameters& args_) = 0;

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

	static CompressionSettings GetCompressionSettings(const InputParameters& args_)
	{
		CompressionSettings settings;

		settings.lossy = args_.lossyCompression;
		settings.dnaOrder = args_.dnaCompressionLevel * 3;

		if (settings.lossy)
			settings.qualityOrder = args_.qualityCompressionLevel * 3;
		else
			settings.qualityOrder = args_.qualityCompressionLevel;

		settings.tagPreserveFlags = args_.tagPreserveFlags;
		settings.calculateCrc32 = args_.calculateCrc32;

		return settings;
	}
};

class DsrcCompressorST : public IDsrcOperator
{
public:
	bool Process(const InputParameters& args_);
};

class DsrcDecompressorST : public IDsrcOperator
{
public:
	bool Process(const InputParameters& args_);
};

class DsrcCompressorMT : public IDsrcOperator
{
public:
	bool Process(const InputParameters& args_);
};

class DsrcDecompressorMT : public IDsrcOperator
{
public:
	bool Process(const InputParameters& args_);
};

} // namespace comp

} // namespace dsrc

#endif
