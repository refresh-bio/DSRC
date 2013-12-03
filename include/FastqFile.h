/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQFILE
#define H_FASTQFILE

#include "Globals.h"
#include "FastqRecord.h"

namespace dsrc
{

namespace wrap
{

class IFastqStream
{
public:
	enum StreamMode
	{
		ModeNone = 0,
		ModeRead,
		ModeWrite
	};

	IFastqStream();
	virtual ~IFastqStream();

	bool ReadNextRecord(FastqRecord& rec_)
	{
		return ReadString(rec_.tag) && ReadString(rec_.sequence)
				&& ReadString(rec_.plus) && ReadString(rec_.quality);
	}

	void WriteNextRecord(const FastqRecord& rec_)
	{
		WriteString(rec_.tag);
		WriteString(rec_.sequence);
		WriteString(rec_.plus);
		WriteString(rec_.quality);
	}

protected:
	struct BufferImpl;
	static const uint64 DefaultBufferSize = 1 << 12;

	StreamMode mode;
	BufferImpl* ioBuffer;

	virtual uint64 ReadBuffer(byte* mem_, uint64 size_) = 0;
	virtual uint64 WriteBuffer(byte* mem_, uint64 size_) = 0;

	void Open(StreamMode mode_);
	void Close();

	bool ReadString(std::string& str_);
	void WriteString(const std::string& str_);
	void Flush();
};

class FastqFile : public IFastqStream
{
public:
	FastqFile();
	~FastqFile();

	void Open(const std::string& filename_);
	void Create(const std::string& filename_);
	void Close();

private:
	struct StreamImpl;

	StreamImpl* impl;

	uint64 ReadBuffer(byte* mem_, uint64 size_);
	uint64 WriteBuffer(byte* mem_, uint64 size_);
};


} // namespace wrap

} // namespace dsrc

#endif // H_FASTQFILE
