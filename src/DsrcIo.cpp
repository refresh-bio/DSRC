/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "DsrcIo.h"
#include "ErrorHandler.h"

namespace dsrc
{

namespace comp
{

void DsrcWriter::operator()()
{
	int64 partId = 0;
	int64 lastPartId = -1;

	DsrcDataChunk* part = NULL;
	std::map<int64, DsrcDataChunk*> partsQueue;

	while (!errorHandler.IsError() && dsrcQueue.Pop(partId, part))
	{
		ASSERT(part->size < (1 << 30));

		if (partId != lastPartId + 1)
		{
			ASSERT(partsQueue.count(partId) == 0);
			partsQueue.insert(std::make_pair(partId, part));
			continue;
		}

		dsrcWriter.WriteNextChunk(part);

		lastPartId++;

		dsrcPool.Release(part);
		part = NULL;

		while (partsQueue.size() > 0)
		{
			partId = partsQueue.begin()->first;
			if (partId == lastPartId + 1)
			{
				part = partsQueue.begin()->second;

				dsrcWriter.WriteNextChunk(part);

				lastPartId++;

				dsrcPool.Release(part);
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
		partId = partsQueue.begin()->first;
		if (partId == lastPartId + 1)
		{
			part = partsQueue.begin()->second;

			dsrcWriter.WriteNextChunk(part);

			lastPartId++;

			dsrcPool.Release(part);
			partsQueue.erase(partsQueue.begin());
		}
		else
		{
			break;
		}
	}

	ASSERT(errorHandler.IsError() || partsQueue.size() == 0);
}

void DsrcReader::operator()()
{
	int64 partId = 0;
	DsrcDataChunk* part = NULL;

	dsrcPool.Acquire(part);

	// waste of last SB --> TODO: fix it, eg. while !EOF
	while (!errorHandler.IsError() && dsrcReader.ReadNextChunk(part))
	{
		ASSERT(part->size < (1 << 30));

		dsrcQueue.Push(partId++, part);
		dsrcPool.Acquire(part);
	}

	dsrcQueue.SetCompleted();
}

} // namespace comp

} // namespace dsrc
