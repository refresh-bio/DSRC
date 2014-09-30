/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQPARSER
#define H_FASTQPARSER

#include "../include/dsrc/Globals.h"

#include "Common.h"
#include "Fastq.h"

namespace dsrc
{

namespace fq
{

class FastqParser
{
public:
	FastqParser();

	uint64 ParseFrom(const FastqDataChunk& chunk_, std::vector<FastqRecord>& records_,
					 uint64& rec_count_, StreamsInfo& streamsInfo_);
	bool Analyze(const FastqDataChunk& chunk_, FastqDatasetType& header_, bool estimateQualityOffset_ = false);

protected:
	core::Buffer* buffer;
	byte* memory;
	uint64 memoryPos;
	uint64 memorySize;
	uint64 skippedBytes;

	bool ReadNextRecord(FastqRecord& rec_)
	{
		if (memoryPos == memorySize)
			return false;

		rec_.title = memory + memoryPos;
		rec_.titleLen = SkipLine();
		if (rec_.titleLen == 0 || rec_.title[0] != '@')
			return false;

		rec_.sequence = memory + memoryPos;
		rec_.sequenceLen = SkipLine();

		// read plus
		uint32 plusLen = SkipLine();

		rec_.quality = memory + memoryPos;
		rec_.qualityLen = SkipLine();

		return (plusLen > 0 && rec_.sequenceLen == rec_.qualityLen);
	}

	bool ReadLine(uchar *str_, uint32& len_, uint32& size_)
	{
		uint32 i = 0;
		for (;;)
		{
			int32 c = Getc();
			if (c == -1)
				break;

			if (c != '\n' && c != '\r')
			{
				if (i >= size_)
				{
					core::extend_string(str_, size_);
				}
				str_[i++] = (uchar)c;
			}
			else
			{
				if (c == '\r' && Peekc() == '\n')	// case of CR LF
					Skipc();

				if (i > 0)
					break;
			}
		}
		str_[i] = 0;
		len_ = i;
		return i > 0;
	}

	uint32 SkipLine()
	{
		uint32 len = 0;
		for (;;)
		{
			int32 c = Getc();
			if (c == -1)
				break;

			if (c != '\n' && c != '\r')
			{
				len++;
			}
			else
			{
				if (c == '\r' && Peekc() == '\n')	// case of CR LF
					Skipc();

				break;
			}
		}
		return len;
	}

	int32 Getc()
	{
		if (memoryPos == memorySize)
			return -1;
		return memory[memoryPos++];
	}

	void Skipc()
	{
		memoryPos++;
		skippedBytes++;
	}

	int32 Peekc()
	{
		if (memoryPos == memorySize)
			return -1;
		return memory[memoryPos];
	}

	void ExtendBuffer(uint32 size_)
	{
		if (buffer->Size() > size_)
			return;

		buffer->Extend(size_, true);
		memory = buffer->Pointer();
		memorySize = buffer->Size();
	}
};


class FastqParserExt : public FastqParser
{
public:
	FastqParserExt()
		:	totalBytesCut(0)
	{}

	// this is a bad design, but trying to avoid virtual function call
	uint64 ParseFrom(const FastqDataChunk &chunk_, std::vector<FastqRecord> &records_, uint64 &rec_count_,
					 StreamsInfo& streamsInfo_, uint64 tagPreserveFlags_);
	bool ReadNextRecord(FastqRecord& rec_, uchar* tagBuffer_, uint64 tagPreserveFlags_);

private:
	static const uint32 MaxTagBufferSize = 512;

	uint64 totalBytesCut;
};

} // namespace fq

} // namespace dsrc


#endif // H_FASTQPARSER
