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
#include "FastqStream.h"


namespace dsrc
{

namespace ext
{

template <class _TFastqStream>
struct TFastqFileImplBase
{
	_TFastqStream* stream;

	TFastqFileImplBase()
		:	stream(NULL)
	{}

	virtual ~TFastqFileImplBase()
	{}

	virtual bool Open(const std::string& filename_)
	{
		if (stream != NULL)
			return false;

		try
		{
			stream = new _TFastqStream(filename_);
		}
		catch (const DsrcException& )
		{
			return false;
		}

		return true;
	}

	virtual void Close()
	{
		if (stream != NULL)
		{
			stream->Close();
			delete stream;
			stream = NULL;
		}
	}
};


// FASTQ file writer/reader implementation
//
FastqFileReader::FastqFileReader()
{}

FastqFileReader::~FastqFileReader()
{}

bool FastqFileReader::AnalyzeFastqDatasetType(FastqDatasetType &type_, byte *buffer_, uint64 bufferSize_)
{
	if (bufferSize_ == 0 || buffer_ == NULL)
		return false;

	fq::FastqDataChunk fqChunk(buffer_, bufferSize_);
	fq::FastqParser parser;

	return parser.Analyze(fqChunk, type_, true);
}


FastqFileWriter::FastqFileWriter()
{}

FastqFileWriter::~FastqFileWriter()
{}


// FASTQ file records writer/reader implementation
//
struct FastqFileRecordsReader::FastqRecordsReaderImpl : public TFastqFileImplBase<fq::FastqFileReader>
{
	static const uint64 DefaultBufferSize = 1 << 13;		// 8K
	fq::FastqDataChunk buffer;
	uint64 bufferPos;

	FastqRecordsReaderImpl(uint64 bufferSize_ = DefaultBufferSize)
		:	buffer(bufferSize_)
		,	bufferPos(0)
	{}

	bool ReadString(std::string& str_)
	{
		str_.clear();

		if (bufferPos >= buffer.size)
		{
			FeedBuffer();
			if (buffer.size == 0)
				return false;
		}

		const char* data = (const char*)buffer.data.Pointer();
		for ( ; ; )
		{
			ASSERT(bufferPos < buffer.size);

			char c = data[bufferPos++];
			if (c == '\n')
				break;

			str_.push_back(c);

			// we are sure that whole blocks of FASTQ file are read
			if (bufferPos >= buffer.size)
				break;
		}

		return str_.length() > 0;
	}

	bool FeedBuffer()
	{
		stream->ReadNextChunk(&buffer);
		bufferPos = 0;

		return buffer.size > 0;
	}
};

struct FastqFileRecordsWriter::FastqRecordsWriterImpl: public TFastqFileImplBase<fq::FastqFileWriter>
{
	static const uint64 DefaultBufferSize = 1 << 13;		// 8K
	fq::FastqDataChunk buffer;
	uint64 bufferPos;

	FastqRecordsWriterImpl(uint64 bufferSize_ = DefaultBufferSize)
		:	buffer(bufferSize_)
		,	bufferPos(0)
	{}

	bool WriteString(const std::string& str_)
	{
		char* data = (char*)buffer.data.Pointer();
		if (str_.length() + bufferPos + 1 > buffer.data.Size())
		{
			// flush all the buffer data
			FlushBuffer();
		}

		std::copy(str_.c_str(), str_.c_str() + str_.length(), data + bufferPos);

		bufferPos += str_.length();
		data[bufferPos++] = '\n';

		buffer.size += str_.length() + 1;

		return true;
	}

	void FlushBuffer()
	{
		if (bufferPos > 0)
		{
			stream->WriteNextChunk(&buffer);	// TODO: implement proper error handling on write
			bufferPos = 0;
			buffer.size = 0;
		}
	}

	void Close()
	{
		FlushBuffer();

		TFastqFileImplBase<fq::FastqFileWriter>::Close();
	}
};


/*
bool FastqFile::GetFastqDatasetType(FastqDatasetType &type_)
{
	if (fastqImpl->mode != FastqFileImpl::ModeRead)
		return false;

	if (fastqImpl->bufferPos == 0 && fastqImpl->bufferSize == 0)
		fastqImpl->FeedBuffer();

	return FastqFile::AnalyzeFastqDatasetType(type_, fastqImpl->buffer.Pointer(), fastqImpl->buffer.Size());
}
*/

// FASTQ file records writer/reader class
//
FastqFileRecordsReader::FastqFileRecordsReader()
{
	readerImpl = new FastqRecordsReaderImpl();
}

FastqFileRecordsReader::~FastqFileRecordsReader()
{
	delete readerImpl;
}

bool FastqFileRecordsReader::Open(const std::string &filename_)
{
	return readerImpl->Open(filename_);
}

void FastqFileRecordsReader::Close()
{
	readerImpl->Close();
}

bool FastqFileRecordsReader::ReadNextRecord(FastqRecord& rec_)
{
	if (readerImpl->stream == NULL)
		return false;

	return readerImpl->ReadString(rec_.tag)
			&& readerImpl->ReadString(rec_.sequence)
			&& readerImpl->ReadString(rec_.plus)
			&& readerImpl->ReadString(rec_.quality);
}


FastqFileRecordsWriter::FastqFileRecordsWriter()
{
	writerImpl = new FastqRecordsWriterImpl();
}

FastqFileRecordsWriter::~FastqFileRecordsWriter()
{
	delete writerImpl;
}

bool FastqFileRecordsWriter::Open(const std::string &filename_)
{
	return writerImpl->Open(filename_);
}

void FastqFileRecordsWriter::Close()
{
	writerImpl->Close();
}

bool FastqFileRecordsWriter::WriteNextRecord(const FastqRecord& rec_)
{
	if (writerImpl->stream == NULL)
		return false;

	return writerImpl->WriteString(rec_.tag)
			&& writerImpl->WriteString(rec_.sequence)
			&& writerImpl->WriteString(rec_.plus)
			&& writerImpl->WriteString(rec_.quality);
}


// FASTQ file records writer/reader implementation -- TODO: bad design, too much code bloat...
//
struct FastqFileBlocksReader::FastqBlocksReaderImpl : public TFastqFileImplBase<fq::FastqFileReader>
{};

struct FastqFileBlocksWriter::FastqBlocksWriterImpl : public TFastqFileImplBase<fq::FastqFileWriter>
{};


// FASTQ file records writer/reader class
//
FastqFileBlocksReader::FastqFileBlocksReader()
{
	readerImpl = new FastqBlocksReaderImpl();
}

FastqFileBlocksReader::~FastqFileBlocksReader()
{
	delete readerImpl;
}

bool FastqFileBlocksReader::Open(const std::string &filename_)
{
	return readerImpl->Open(filename_);
}

void FastqFileBlocksReader::Close()
{
	readerImpl->Close();
}

unsigned long int FastqFileBlocksReader::ReadNextBlock(char* buffer_, unsigned long int bufferSize_)
{
	if (readerImpl->stream == NULL || buffer_ == NULL || bufferSize_ <= MinBufferSize)
		return 0;

	fq::FastqDataChunk chunk((byte*)buffer_, (uint64)bufferSize_);

	readerImpl->stream->ReadNextChunk(&chunk);
	return (unsigned long int)chunk.size;
}


FastqFileBlocksWriter::FastqFileBlocksWriter()
{
	writerImpl = new FastqBlocksWriterImpl();
}

FastqFileBlocksWriter::~FastqFileBlocksWriter()
{
	delete writerImpl;
}

bool FastqFileBlocksWriter::Open(const std::string &filename_)
{
	return writerImpl->Open(filename_);
}

void FastqFileBlocksWriter::Close()
{
	writerImpl->Close();
}

unsigned long int FastqFileBlocksWriter::WriteNextBlock(char* buffer_, unsigned long int dataSize_)
{
	if (writerImpl->stream == NULL || buffer_ == NULL || dataSize_ == 0)
		return 0;

	fq::FastqDataChunk chunk((byte*)buffer_, (uint64)dataSize_);
	chunk.size = dataSize_;

	writerImpl->stream->WriteNextChunk(&chunk);		// TODO: implement proper error handling on write
	return dataSize_;
}


} // namespace ext

} // namespace dsrc
