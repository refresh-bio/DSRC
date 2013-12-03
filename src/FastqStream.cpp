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
	const uint64 cbuf_size = chunk_->data.Size();
	chunk_->size = 0;

	int64 to_read = cbuf_size - bufferSize;
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = Read(data + chunk_->size, to_read);

	if (r > 0)
	{
		if (r == to_read)	// somewhere before end
		{
			uint64 chunk_end = cbuf_size - SwapBufferSize;

			chunk_end = GetNextRecordPos(data, chunk_end, cbuf_size);

			chunk_->size = chunk_end - 1;

			std::copy(data + chunk_end, data + cbuf_size, swapBuffer.Pointer());
			bufferSize = cbuf_size - chunk_end;
		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
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

