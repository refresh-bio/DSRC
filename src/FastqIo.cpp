/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include "FastqIo.h"

#include <vector>
#include <map>

#include "Buffer.h"
#include "FastqStream.h"
#include "FastqParser.h"
#include "ErrorHandler.h"

namespace dsrc
{

namespace fq
{

bool FastqReader::AnalyzeFirstChunk(FastqDatasetType& header_, bool estimateQualityOffset_)
{
	FastqParser parser;
	FastqDataChunk* fqChunk;
	recordsPool.Acquire(fqChunk);

	// analyze first chunk
	if (!fileReader.ReadNextChunk(fqChunk) || !parser.Analyze(*fqChunk, header_, estimateQualityOffset_))
	{
		recordsPool.Release(fqChunk);
		return false;
	}

	// add the chunk to later process queue
	ASSERT(numParts == 0);
	recordsQueue.Push(numParts++, fqChunk);

	return true;
}

void FastqReader::operator()()
{
	FastqDataChunk* part = NULL;

	recordsPool.Acquire(part);

	while (!errorHandler.IsError() && fileReader.ReadNextChunk(part))
	{
		ASSERT(part->size > 0);

		recordsQueue.Push(numParts, part);
		numParts++;

		recordsPool.Acquire(part);
	}

	ASSERT(part->size == 0);
	recordsPool.Release(part);		// the last empty part

	recordsQueue.SetCompleted();
}


// FastqWriter
//
void FastqWriter::operator()()
{
	FastqDataChunk* part = NULL;
	int64 partId = 0;
	int64 lastPartId = -1;

	std::map<int64, FastqDataChunk*> partsQueue;

	while (!errorHandler.IsError() && recordsQueue.Pop(partId, part))
	{
		ASSERT(part->size > 0);
		std::vector<FastqRecord>::const_iterator i;

		if (partId != lastPartId + 1)
		{
			ASSERT(partsQueue.count(partId) == 0);
			partsQueue.insert(std::make_pair(partId, part));
			continue;
		}

		fileWriter.WriteNextChunk(part);
		lastPartId++;

		recordsPool.Release(part);
		part = NULL;

		while (partsQueue.size() > 0)
		{
			if (partsQueue.begin()->first == lastPartId + 1)
			{
				part = partsQueue.begin()->second;
				fileWriter.WriteNextChunk(part);
				lastPartId++;

				recordsPool.Release(part);
				partsQueue.erase(partsQueue.begin());
			}
			else
			{
				break;
			}
		}

		part = NULL;
	}

	while (!errorHandler.IsError() && partsQueue.size() > 0)
	{
		if (partsQueue.begin()->first == lastPartId + 1)
		{
			part = partsQueue.begin()->second;
			fileWriter.WriteNextChunk(part);
			lastPartId++;

			recordsPool.Release(part);
			partsQueue.erase(partsQueue.begin());
		}
		else
		{
			break;
		}
	}

	ASSERT(errorHandler.IsError() || partsQueue.size() == 0);
}

} // namesapce fq

} // namespace dsrc
