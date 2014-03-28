/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "FastqStream.h"

namespace dsrc
{

namespace fq
{

bool IFastqStreamReader::ReadNextChunk(FastqDataChunk* chunk_)
{
	if (Eof())
	{
		chunk_->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	uchar* data = chunk_->data.Pointer();
	const uint64 cbufSize = chunk_->data.Size();
	chunk_->size = 0;

	int64 toRead = cbufSize - bufferSize;
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = Read(data + chunk_->size, toRead);

	if (r > 0)
	{
		if (r == toRead)	// somewhere before end
		{
			uint64 chunkEnd = cbufSize - SwapBufferSize;

			chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);

			chunk_->size = chunkEnd - 1;
			if (usesCrlf)
				chunk_->size -= 1;

			std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
			bufferSize = cbufSize - chunkEnd;
		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
			if (usesCrlf)
				chunk_->size -= 1;

			eof = true;
		}
	}
	else
	{
		eof = true;
	}

	return true;
}

uint64 IFastqStreamReader::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
{
	SkipToEol(data_, pos_, size_);
	++pos_;

	// find beginning of the next record
	while (data_[pos_] != '@')
	{
		SkipToEol(data_, pos_, size_);
		++pos_;
	}
	uint64 pos0 = pos_;

	SkipToEol(data_, pos_, size_);
	++pos_;

	if (data_[pos_] == '@')			// previous one was a quality field
		return pos_;

	SkipToEol(data_, pos_, size_);
	++pos_;

	ASSERT(data_[pos_] == '+');	// pos0 was the start of tag
	return pos0;
}

} // namespace fq

} // namespace dsrc

