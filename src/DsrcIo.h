/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCWRITER
#define H_DSRCWRITER

#include "../include/dsrc/Globals.h"

#include "Buffer.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "DsrcFile.h"

#include <vector>
#include <map>

namespace dsrc
{

namespace comp
{

struct DsrcDataChunk : public core::DataChunk
{
	fq::StreamsInfo rawStreamsInfo;
	fq::StreamsInfo compStreamsInfo;

	DsrcDataChunk(uint64 bufferSize_ = core::DataChunk::DefaultBufferSize)
		:	core::DataChunk(bufferSize_)
	{}

	void Reset()
	{
		core::DataChunk::Reset();
		rawStreamsInfo.Clear();
		compStreamsInfo.Clear();
	}
};

typedef core::TDataQueue<DsrcDataChunk> DsrcDataQueue;
typedef core::TDataPool<DsrcDataChunk> DsrcDataPool;

class IDsrcIoOperator
{
public:
	IDsrcIoOperator(DsrcDataQueue& queue_, DsrcDataPool& pool_, core::ErrorHandler& errorHandler_)
		:	dsrcQueue(queue_)
		,	dsrcPool(pool_)
		,	errorHandler(errorHandler_)
	{}

	virtual ~IDsrcIoOperator() {}

	virtual void operator()() = 0;

protected:
	DsrcDataQueue& dsrcQueue;
	DsrcDataPool& dsrcPool;
	core::ErrorHandler& errorHandler;
};

class DsrcWriter : public IDsrcIoOperator
{
	DsrcFileWriter& dsrcWriter;

public:
	DsrcWriter(DsrcFileWriter& writer_, DsrcDataQueue& queue_, DsrcDataPool& pool_, core::ErrorHandler& errorHandler_)
		:	IDsrcIoOperator(queue_, pool_, errorHandler_)
		,	dsrcWriter(writer_)
	{}

	void operator()();
};


class DsrcReader : public IDsrcIoOperator
{
	DsrcFileReader& dsrcReader;

public:
	DsrcReader(DsrcFileReader& reader_, DsrcDataQueue& queue_, DsrcDataPool& pool_, core::ErrorHandler& errorHandler_)
		:	IDsrcIoOperator(queue_, pool_, errorHandler_)
		,	dsrcReader(reader_)
	{}

	void operator()();
};

} // namespace comp

} // namespace dsrc

#endif
