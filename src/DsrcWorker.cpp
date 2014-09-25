/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include "../include/dsrc/Globals.h"

#include "DsrcWorker.h"
#include "DsrcIo.h"
#include "BlockCompressor.h"
#include "ErrorHandler.h"

#include <algorithm>
#include <cstring>


namespace dsrc
{

namespace comp
{

using namespace core;
using namespace fq;

void DsrcCompressor::Process()
{
	int64 partId = 0;

	FastqDataChunk* fqChunk = NULL;
	DsrcDataChunk* dsrcData = NULL;

	BlockCompressor superblock(datasetType, compSettings);

	while (!errorHandler.IsError() && fastqQueue.Pop(partId, fqChunk))
	{
		ASSERT(fqChunk->size > 0);

		dsrcPool.Acquire(dsrcData);
		ASSERT(dsrcData != NULL);

		BitMemoryWriter bitMemory(dsrcData->data);

		superblock.Store(bitMemory, dsrcData->rawStreamsInfo, dsrcData->compStreamsInfo, *fqChunk);

		bitMemory.Flush();
		dsrcData->size = bitMemory.Position();

		if (compSettings.calculateCrc32)
		{
			BitMemoryReader reader(dsrcData->data.Pointer(), dsrcData->data.Size());
			std::fill(fqChunk->data.Pointer(), fqChunk->data.Pointer() + fqChunk->data.Size(), 0xCC);

			if (!superblock.VerifyChecksum(reader, *fqChunk))
			{
				errorHandler.SetError("CRC32 checksums mismatch.");
			}
		}

		dsrcQueue.Push(partId, dsrcData);
		dsrcData = NULL;
		bitMemory.Reset();

		fastqPool.Release(fqChunk);
		fqChunk = NULL;
	}

	dsrcQueue.SetCompleted();
}

void DsrcDecompressor::Process()
{
	int64 partId = 0;

	FastqDataChunk* fqChunk = NULL;
	DsrcDataChunk* dsrcData = NULL;

	BlockCompressor superblock(datasetType, compSettings);

	while (!errorHandler.IsError() && dsrcQueue.Pop(partId, dsrcData))
	{
		ASSERT(dsrcData);
		ASSERT(dsrcData->size > 0);
		ASSERT(dsrcData->size <= dsrcData->data.Size());

		BitMemoryReader bitMemory(dsrcData->data.Pointer(), dsrcData->size);

		fastqPool.Acquire(fqChunk);

		superblock.Read(bitMemory, *fqChunk);

		fastqQueue.Push(partId, fqChunk);
		fqChunk = NULL;

		dsrcPool.Release(dsrcData);
		dsrcData = NULL;
	}

	fastqQueue.SetCompleted();
}

} // namespace comp

} // namespace dsrc
