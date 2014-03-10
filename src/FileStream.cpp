/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "FileStream.h"

#if defined (_WIN32)
#	define _CRT_SECURE_NO_WARNINGS
#	pragma warning(disable : 4996) // D_SCL_SECURE
#	pragma warning(disable : 4244) // conversion uint64 to uint32
//#	pragma warning(disable : 4267)
#	define FOPEN	fopen
#	define FSEEK	_fseeki64
#	define FTELL	_ftelli64
#	define FCLOSE	fclose
#elif __APPLE__	// Apple by default suport 64 bit file operations (Darwin 10.5+)
#	define FOPEN	fopen
#	define FSEEK	fseek
#	define FTELL	ftell
#	define FCLOSE	fclose
#else
#	if !defined(_LARGEFILE_SOURCE)
#		define _LARGEFILE_SOURCE
#		if !defined(_LARGEFILE64_SOURCE)
#			define _LARGEFILE64_SOURCE
#		endif
#	endif
#	if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#		undef _FILE_OFFSET_BITS
#	endif
#	if !defined(_FILE_OFFSET_BITS)
#		define _FILE_OFFSET_BITS 64
#	endif
#	define FOPEN	fopen64
#	define FSEEK	fseeko64
#	define FTELL	ftello64
#	define FCLOSE	fclose
#endif

#include <stdio.h>

namespace dsrc
{

namespace core
{

struct IFileStream::FileStreamImpl
{
	FILE* file;

	FileStreamImpl()
		:	file(NULL)
	{}
};

IFileStream::IFileStream()
{
	impl = new FileStreamImpl();
}

IFileStream::~IFileStream()
{
	ASSERT(impl != NULL);
	delete impl;
}

FileStreamReader::FileStreamReader(const std::string& fileName_)
{
	FILE* f = FOPEN(fileName_.c_str(), "rb");
	if (f == NULL)
	{
		throw DsrcException(("Cannot open file to read:" + fileName_).c_str());
	}
	impl->file = f;
}

FileStreamReader::~FileStreamReader()
{
	ASSERT(impl != NULL);

	if (impl->file != NULL)
		FCLOSE(impl->file);
}

void FileStreamReader::Close()
{
	ASSERT(impl->file != NULL);

	FCLOSE(impl->file);
	impl->file = NULL;
}

int64 FileStreamReader::Read(uchar *mem_, uint64 size_)
{
	int64 n = fread(mem_, 1, size_, impl->file);
	return n;
}

FileStreamReaderExt::FileStreamReaderExt(const std::string& fileName_)
	:	FileStreamReader(fileName_)
	,	size(0)
	,	position(0)
{
	FSEEK(impl->file, 0, SEEK_END);
	size = FTELL(impl->file);

	FSEEK(impl->file, 0, SEEK_SET);
	position = 0;
}

int64 FileStreamReaderExt::Read(uchar *mem_, uint64 size_)
{
	int64 n = fread(mem_, 1, size_, impl->file);
	if (n >= 0)
		position += n;
	return n;
}

void FileStreamReaderExt::SetPosition(uint64 pos_)
{
	ASSERT(impl->file != NULL);

	if (pos_ > size)
	{
		throw DsrcException("Position exceeds stream size");
	}
	FSEEK(impl->file, pos_, SEEK_SET);
	position = pos_;
}

FileStreamWriter::FileStreamWriter(const std::string& fileName_)
{
	FILE* f = FOPEN(fileName_.c_str(), "wb");
	if (f == NULL)
	{
		throw DsrcException(("Cannot open file to write:" + fileName_).c_str());
	}
	impl->file = f;
}

FileStreamWriter::~FileStreamWriter()
{
	ASSERT(impl != NULL);

	if (impl->file != NULL)
		FCLOSE(impl->file);
}

void FileStreamWriter::Close()
{
	ASSERT(impl->file != NULL);

	FCLOSE(impl->file);
	impl->file = NULL;
}

int64 FileStreamWriter::Write(const uchar *mem_, uint64 size_)
{
	int64 n = fwrite(mem_, 1, size_, impl->file);
	return n;
}

FileStreamWriterExt::FileStreamWriterExt(const std::string& fileName_)
	:	FileStreamWriter(fileName_)
	,	position(0)
{}

int64 FileStreamWriterExt::Write(const uchar *mem_, uint64 size_)
{
	int64 n = fwrite(mem_, 1, size_, impl->file);
	if (n >= 0)
		position += n;
	return n;
}

void FileStreamWriterExt::SetPosition(uint64 pos_)
{
	ASSERT(impl->file != NULL);

	FSEEK(impl->file, pos_, SEEK_SET);
	position = pos_;
}

} // namespace core

} // namespace dsrc
