/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef CPP_CONFIGURABLE
#define CPP_CONFIGURABLE

#include "../include/dsrc/Configurable.h"

#include "Common.h"

namespace dsrc
{

namespace wrap
{

struct Configurable::ConfigImpl
{
	comp::InputParameters inputParams;
	bool plusRepetition;	// optional
	bool colorSpace;		// optional

	ConfigImpl()
		:	plusRepetition(false)
		,	colorSpace(false)
	{}

	static ConfigImpl Default()
	{
		ConfigImpl pars;
		pars.inputParams = comp::InputParameters::Default();
		pars.plusRepetition = false;
		pars.colorSpace = false;
		return pars;
	}
};


Configurable::Configurable()
	:	config(NULL)
{
	config = new ConfigImpl;
}

Configurable::~Configurable()
{
	delete config;
}

void Configurable::SetFastqBufferSizeMB(uint64 size_)
{
	if (size_ == 0 || size_ > 1024)
		throw DsrcException("Invalid argument: invalid FASTQ buffer size [1-1024]");

	config->inputParams.fastqBufferSizeMB = size_;
}

uint64 Configurable::GetFastqBufferSizeMB() const
{
	return config->inputParams.fastqBufferSizeMB;
}

void Configurable::SetDnaCompressionLevel(uint32 level_)
{
	if (level_ > 3)
		throw DsrcException("Invalid argument: invalid DNA compression level [0-3]");

	config->inputParams.dnaCompressionLevel = level_;
}

uint32 Configurable::GetDnaCompressionLevel() const
{
	return config->inputParams.dnaCompressionLevel;
}

void Configurable::SetQualityCompressionLevel(uint32 level_)
{
	if (level_ > 2)
		throw DsrcException("Invalid argument: invalid Quality compression level [0-2]");

	config->inputParams.qualityCompressionLevel = level_;
}

uint32 Configurable::GetQualityCompressionLevel() const
{
	return config->inputParams.qualityCompressionLevel;
}

void Configurable::SetLossyCompression(bool lossy_)
{
	config->inputParams.lossyCompression = lossy_;
}

bool Configurable::IsLossyCompression() const
{
	return config->inputParams.lossyCompression;
}

void Configurable::SetQualityOffset(uint32 off_)
{
	if (off_ != 33 && off_ != 64)
		throw DsrcException("Invalid argument: only valid Quality offset are 33 and 64");

	config->inputParams.qualityOffset = off_;
}

uint32 Configurable::GetQualityOffset() const
{
	return config->inputParams.qualityOffset;
}

void Configurable::SetColorSpace(bool cs_)
{
	config->colorSpace = cs_;
}

bool Configurable::IsColorSpace() const
{
	return config->colorSpace;
}

void Configurable::SetPlusRepetition(bool rep_)
{
	config->plusRepetition = rep_;
}

bool Configurable::IsPlusRepetition() const
{
	return config->plusRepetition;
}

void Configurable::SetThreadsNumber(uint32 threadNum_)
{
	if (threadNum_ == 0)
		throw DsrcException("Invalid argument: thread number must be greater than 0");

	config->inputParams.threadNum = threadNum_;
}

uint32 Configurable::GetThreadsNumber() const
{
	return config->inputParams.threadNum;
}

void Configurable::SetStdIoUsing(bool use_)
{
	config->inputParams.useFastqStdIo = use_;
}

bool Configurable::IsStdIoUsing() const
{
	return config->inputParams.useFastqStdIo;
}

void Configurable::SetCrc32Checking(bool use_)
{
	config->inputParams.calculateCrc32 = use_;
}

bool Configurable::IsCrc32Checking() const
{
	return config->inputParams.calculateCrc32;
}

void Configurable::SetTagFieldFilterMask(uint64 mask_)
{
	config->inputParams.tagPreserveFlags = mask_;
}

uint64 Configurable::GetTagFieldFilterMask() const
{
	return config->inputParams.tagPreserveFlags;
}

void* Configurable::GetInputParameters() const
{
	return (void*)&config->inputParams;
}

} // namespace wrap

} // namespace dsrc


#endif // CPP_CONFIGURABLE
