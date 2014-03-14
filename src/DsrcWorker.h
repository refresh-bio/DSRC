/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCWORKER
#define H_DSRCWORKER

#include "../include/dsrc/Globals.h"

#include "Fastq.h"
#include "BitMemory.h"
#include "FastqIo.h"
#include "DsrcIo.h"

#include <vector>
#include <string>

namespace dsrc
{

namespace comp
{

class IDsrcThreadWorker
{
public:
	IDsrcThreadWorker(fq::FastqDataQueue& fastqQueue_, fq::FastqDataPool& fastqPool_,
				   DsrcDataQueue& dsrcQueue_, DsrcDataPool& dsrcPool_, core::ErrorHandler& errorHandler_,
				   const fq::FastqDatasetType& type_, const CompressionSettings& settings_)
		:	fastqQueue(fastqQueue_)
		,	fastqPool(fastqPool_)
		,	dsrcQueue(dsrcQueue_)
		,	dsrcPool(dsrcPool_)
		,	errorHandler(errorHandler_)
		,	datasetType(type_)
		,	compSettings(settings_)
	{}

	virtual ~IDsrcThreadWorker() {}

	void operator() ()
	{
		Process();
	}

protected:
	fq::FastqDataQueue&	fastqQueue;
	fq::FastqDataPool& fastqPool;
	DsrcDataQueue& dsrcQueue;
	DsrcDataPool& dsrcPool;
	core::ErrorHandler&	errorHandler;

	fq::FastqDatasetType datasetType;
	CompressionSettings compSettings;

private:
	virtual void Process() = 0;
};

class DsrcCompressor : public IDsrcThreadWorker
{
public:
	DsrcCompressor(fq::FastqDataQueue& fastqQueue_, fq::FastqDataPool& fastqPool_,
				   DsrcDataQueue& dsrcQueue_, DsrcDataPool& dsrcPool_, core::ErrorHandler& errorHandler_,
				   const fq::FastqDatasetType& type_, const CompressionSettings& settings_)
		:	IDsrcThreadWorker(fastqQueue_, fastqPool_, dsrcQueue_, dsrcPool_, errorHandler_, type_, settings_)
	{}

private:
	void Process();
};


class DsrcDecompressor : public IDsrcThreadWorker
{
public:
	DsrcDecompressor(fq::FastqDataQueue& fastqQueue_, fq::FastqDataPool& fastqPool_,
					DsrcDataQueue& dsrcQueue_, DsrcDataPool& dsrcPool_, core::ErrorHandler& errorHandler_,
					const fq::FastqDatasetType& type_, const CompressionSettings& settings_)
		:	IDsrcThreadWorker(fastqQueue_, fastqPool_, dsrcQueue_, dsrcPool_, errorHandler_, type_, settings_)
	{}

private:
	void Process();
};

} // namespace comp

} // namespace dsrc

#endif
