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
#include "FastqParser.h"


namespace dsrc
{

namespace ext
{

struct FastqFile::FastqFileImpl
{
	static const uint64 DefaultBufferSize = 1 << 12;

	enum StreamMode
	{
		ModeNone = 0,
		ModeRead,
		ModeWrite
	};

	core::IDataStream* stream;
	core::Buffer buffer;
	uint64 size;
	uint64 pos;
	StreamMode mode;

	FastqFileImpl(uint64 bufferSize_ = DefaultBufferSize)
		:	stream(NULL)
		,	buffer(bufferSize_)
		,	size(0)
		,	pos(0)
		,	mode(ModeNone)
	{}

	bool Open(StreamMode mode_)
	{
		if (mode == ModeNone)
			return false;

		mode = mode_;
		size = 0;
		pos = 0;
		return true;
	}

	void Close()
	{
		mode = ModeNone;
		size = 0;
		pos = 0;
	}

	uint64 ReadBuffer(byte* mem_, uint64 size_)
	{
		ASSERT(mode == ModeRead);
		return stream->PerformIo(mem_, size_);
	}

	uint64 WriteBuffer(byte* mem_, uint64 size_)
	{
		ASSERT(mode == ModeWrite);
		return stream->PerformIo(mem_, size_);
	}

	bool ReadString(std::string& str_)
	{
		str_.clear();

		const char* data = (const char*)buffer.Pointer();
		for ( ; ; )
		{
			if (pos >= size)
			{
				FeedBuffer();

				if (size == 0)
					break;
			}

			ASSERT(pos < buffer.Size());

			char c = data[pos++];
			if (c == '\n')
				break;

			str_.push_back(c);
		}

		return str_.length() > 0;
	}

	void WriteString(const std::string& str_)
	{
		char* data = (char*)buffer.Pointer();
		if (str_.length() + pos + 1 > buffer.Size())
		{
			WriteBuffer(buffer.Pointer(), pos);
			pos = 0;
		}

		std::copy(str_.c_str(), str_.c_str() + str_.length(), data + pos);
		pos += str_.length();
		data[pos++] = '\n';
	}


	void FlushBuffer()
	{
		if (mode == ModeWrite && pos > 0)
		{
			WriteBuffer(buffer.Pointer(), pos);
			pos = 0;
		}
	}

	bool FeedBuffer()
	{
		size = ReadBuffer(buffer.Pointer(), buffer.Size());
		pos = 0;

		return size > 0;
	}
};


FastqFile::FastqFile()
	:	fastqImpl(NULL)
{
	fastqImpl = new FastqFileImpl();
}

FastqFile::~FastqFile()
{
	delete fastqImpl;
}

bool FastqFile::Open(const std::string &filename_)
{
	if (!fastqImpl->Open(FastqFileImpl::ModeRead))
		return false;

	ASSERT(fastqImpl->stream != NULL);
	fastqImpl->stream = new core::FileStreamReader(filename_.c_str());

	return true;
}

bool FastqFile::Create(const std::string &filename_)
{
	if (!fastqImpl->Open(FastqFileImpl::ModeWrite))
		return false;

	ASSERT(fastqImpl->stream != NULL);
	fastqImpl->stream = new core::FileStreamWriter(filename_.c_str());

	return true;
}

void FastqFile::Close()
{
	if (fastqImpl->mode == FastqFileImpl::ModeNone)
		return;

	ASSERT(fastqImpl->stream != NULL);

	fastqImpl->FlushBuffer();

	fastqImpl->stream->Close();
	delete fastqImpl->stream;
	fastqImpl->stream = NULL;

	fastqImpl->Close();
}

bool FastqFile::ReadNextRecord(FastqRecord& rec_)
{
	if (fastqImpl->mode == FastqFileImpl::ModeNone)
		return false;

	return fastqImpl->ReadString(rec_.tag)
			&& fastqImpl->ReadString(rec_.sequence)
			&& fastqImpl->ReadString(rec_.plus)
			&& fastqImpl->ReadString(rec_.quality);
}

void FastqFile::WriteNextRecord(const FastqRecord& rec_)
{
	if (fastqImpl->mode == FastqFileImpl::ModeNone)
		return;

	fastqImpl->WriteString(rec_.tag);
	fastqImpl->WriteString(rec_.sequence);
	fastqImpl->WriteString(rec_.plus);
	fastqImpl->WriteString(rec_.quality);
}

bool FastqFile::GetFastqDatasetType(FastqDatasetType &type_)
{
	if (fastqImpl->mode != FastqFileImpl::ModeRead)
		return false;

	if (fastqImpl->pos == 0 && fastqImpl->size == 0)
		fastqImpl->FeedBuffer();

	return FastqFile::AnalyzeFastqDatasetType(type_, fastqImpl->buffer.Pointer(), fastqImpl->buffer.Size());
}

bool FastqFile::AnalyzeFastqDatasetType(FastqDatasetType &type_, byte *buffer_, uint64 bufferSize_)
{
	if (bufferSize_ == 0 || buffer_ == NULL)
		return false;

	fq::FastqDataChunk fqChunk(buffer_, bufferSize_);
	fq::FastqParser parser;

	return parser.Analyze(fqChunk, type_, true);
}

} // namespace ext

} // namespace dsrc
