/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "../include/dsrc/FastqFile.h"

#include "Common.h"
#include "Buffer.h"
#include "DataStream.h"
#include "FastqStream.h"

namespace dsrc
{

namespace wrap
{

struct IFastqStream::BufferImpl
{
	core::Buffer buffer;
	uint64 size;
	uint64 pos;

	BufferImpl(uint64 bufferSize_)
		:	buffer(bufferSize_)
		,	size(0)
		,	pos(0)
	{}
};

IFastqStream::IFastqStream()
	:	mode(ModeNone)
	,	ioBuffer(NULL)
{
	ioBuffer = new BufferImpl(DefaultBufferSize);
}

IFastqStream::~IFastqStream()
{
	delete ioBuffer;
}

void IFastqStream::Open(StreamMode mode_)
{
	if (mode != ModeNone)
		throw DsrcException("Invalid state");

	mode = mode_;
	ioBuffer->size = 0;
	ioBuffer->pos = 0;
}

void IFastqStream::Close()
{
	if (mode == ModeNone)
		throw DsrcException("Invalid state");

	mode = ModeNone;
	ioBuffer->size = 0;
	ioBuffer->pos = 0;
}

bool IFastqStream::ReadString(std::string& str_)
{
	str_.clear();

	const char* data = (const char*)ioBuffer->buffer.Pointer();
	for ( ; ; )
	{
		if (ioBuffer->pos >= ioBuffer->size)
		{
			ioBuffer->pos = 0;
			ioBuffer->size = ReadBuffer(ioBuffer->buffer.Pointer(), ioBuffer->buffer.Size());
			ASSERT(ioBuffer->size <= ioBuffer->buffer.Size());

			if (ioBuffer->size == 0)
				break;
		}

		ASSERT(ioBuffer->pos < ioBuffer->buffer.Size());

		char c = data[ioBuffer->pos++];
		if (c == '\n')
			break;

		str_.push_back(c);
	}

	return str_.length() > 0;
}

void IFastqStream::WriteString(const std::string& str_)
{
	char* data = (char*)ioBuffer->buffer.Pointer();
	if (str_.length() + ioBuffer->pos + 1 > ioBuffer->buffer.Size())
	{
		WriteBuffer(ioBuffer->buffer.Pointer(), ioBuffer->pos);
		ioBuffer->pos = 0;
	}

	std::copy(str_.c_str(), str_.c_str() + str_.length(), data + ioBuffer->pos);
	ioBuffer->pos += str_.length();
	data[ioBuffer->pos++] = '\n';
}

void IFastqStream::Flush()
{
	if (mode == ModeWrite && ioBuffer->pos > 0)
	{
		WriteBuffer(ioBuffer->buffer.Pointer(), ioBuffer->pos);
		ioBuffer->pos = 0;
	}
}


struct FastqFile::StreamImpl
{
	core::IDataStream* stream;
};

FastqFile::FastqFile()
	:	IFastqStream()
	,	impl(NULL)
{
	impl = new StreamImpl();
}

FastqFile::~FastqFile()
{
	if (impl->stream != NULL)
		delete impl->stream;
	delete impl;
}

void FastqFile::Open(const std::string &filename_)
{
	IFastqStream::Open(ModeRead);

	if (impl->stream != NULL)
		throw DsrcException("Invalid state");

	impl->stream = new core::FileStreamReader(filename_.c_str());
}

void FastqFile::Create(const std::string &filename_)
{
	IFastqStream::Open(ModeWrite);

	if (impl->stream != NULL)
		throw DsrcException("Invalid state");

	impl->stream = new core::FileStreamWriter(filename_.c_str());
}

void FastqFile::Close()
{
	if (impl->stream == NULL)
		throw DsrcException("Invalid state");

	Flush();

	impl->stream->Close();

	delete impl->stream;
	impl->stream = NULL;

	IFastqStream::Close();
}

uint64 FastqFile::ReadBuffer(byte* mem_, uint64 size_)
{
	ASSERT(mode == ModeRead);

	return impl->stream->PerformIo(mem_, size_);
}

uint64 FastqFile::WriteBuffer(byte* mem_, uint64 size_)
{
	ASSERT(mode == ModeWrite);

	return impl->stream->PerformIo(mem_, size_);
}

} // namespace wrap

} // namespace dsrc
