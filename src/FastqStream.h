/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
#ifndef H_FASTQSTREAM
#define H_FASTQSTREAM

#include "../include/dsrc/Globals.h"

//#include "Common.h"
#include "Fastq.h"
#include "Buffer.h"
#include "FileStream.h"
#include "StdStream.h"


namespace dsrc
{

namespace fq
{

class IFastqStreamReader
{
private:
	static const uint32 SwapBufferSize = 1 << 13;

public:
	IFastqStreamReader()
		:	stream(NULL)
		,	swapBuffer(SwapBufferSize)
		,	bufferSize(0)
		,	eof(false)
		,	usesCrlf(false)
	{}

	virtual ~IFastqStreamReader()
	{}

	bool Eof() const
	{
		return eof;
	}

	bool ReadNextChunk(FastqDataChunk* chunk_);

	void Close()
	{
		ASSERT(stream != NULL);
		stream->Close();
	}

protected:
	core::IDataStreamReader* stream;

	int64 Read(byte* memory_, uint64 size_)
	{
		ASSERT(stream != NULL);
		return stream->Read(memory_, size_);
	}

private:
	core::Buffer	swapBuffer;
	uint64			bufferSize;
	bool			eof;
	bool			usesCrlf;

	uint64 GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_);

	void SkipToEol(uchar* data_, uint64& pos_, const uint64 size_)
	{
		ASSERT(pos_ < size_);

		while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
			++pos_;

		if (data_[pos_] == '\r' && pos_ < size_)
		{
			if (data_[pos_ + 1] == '\n')
			{
				usesCrlf = true;
				++pos_;
			}
		}
	}

};

class IFastqStreamWriter
{
public:
	IFastqStreamWriter()
	{}

	virtual ~IFastqStreamWriter()
	{}

	void WriteNextChunk(const FastqDataChunk* chunk_)
	{
		ASSERT(stream != NULL);
		stream->Write(chunk_->data.Pointer(), chunk_->size);
	}

	void Close()
	{
		ASSERT(stream != NULL);
		stream->Close();
	}

protected:
	core::IDataStreamWriter* stream;
};


// wrappers
//
class FastqFileReader : public IFastqStreamReader
{
public:
	FastqFileReader(const std::string& fileName_)
	{
		stream = new core::FileStreamReader(fileName_);
	}

	~FastqFileReader()
	{
		delete stream;
	}
};

class FastqFileWriter : public IFastqStreamWriter
{
public:
	FastqFileWriter(const std::string& fileName_)
	{
		stream = new core::FileStreamWriter(fileName_);
	}

	~FastqFileWriter()
	{
		delete stream;
	}
};

class FastqStdIoReader : public IFastqStreamReader
{
public:
	FastqStdIoReader()
	{
		stream = new core::StdStreamReader();
	}

	~FastqStdIoReader()
	{
		delete stream;
	}
};

class FastqStdIoWriter : public IFastqStreamWriter
{
public:
	FastqStdIoWriter()
	{
		stream = new core::StdStreamWriter();
	}

	~FastqStdIoWriter()
	{
		delete stream;
	}
};

} // namespace fq

} // namespace dsrc

#endif // H_FASTQSTREAM
