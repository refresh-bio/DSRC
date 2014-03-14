/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQREADER
#define H_FASTQREADER

#include "../include/dsrc/Globals.h"

#include "Common.h"
#include "Fastq.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "FastqStream.h"

namespace dsrc
{

namespace fq
{

typedef core::TDataQueue<FastqDataChunk> FastqDataQueue;
typedef core::TDataPool<FastqDataChunk> FastqDataPool;


// refactor to one interface
//
class IFastqIoOperator
{
public:
	IFastqIoOperator(FastqDataQueue& queue_, FastqDataPool& pool_, core::ErrorHandler& errorHandler_)
		:	recordsQueue(queue_)
		,	recordsPool(pool_)
		,	errorHandler(errorHandler_)
	{}

	virtual ~IFastqIoOperator() {}

	virtual void operator()() = 0;

protected:
	FastqDataQueue&		recordsQueue;
	FastqDataPool&		recordsPool;
	core::ErrorHandler& errorHandler;
};

class FastqReader : public IFastqIoOperator
{
public:
	FastqReader(IFastqStreamReader& reader_, FastqDataQueue& queue_, FastqDataPool& pool_, core::ErrorHandler& errorHandler_)
		:	IFastqIoOperator(queue_, pool_, errorHandler_)
		,	fileReader(reader_)
		,	numParts(0)
	{}

	bool AnalyzeFirstChunk(FastqDatasetType& header_, bool estimateQualityOffset_);
	void operator()();

private:
	IFastqStreamReader&	fileReader;
	uint32 numParts;
};

class FastqWriter : public IFastqIoOperator
{
public:
	FastqWriter(IFastqStreamWriter& writer_, FastqDataQueue& queue_, FastqDataPool& pool_, core::ErrorHandler& errorHandler_)
		:	IFastqIoOperator(queue_, pool_, errorHandler_)
		,	fileWriter(writer_)
	{}

	void operator()();

private:
	IFastqStreamWriter&	fileWriter;
};

} // namespace fq

} // namespace dsrc

#endif
