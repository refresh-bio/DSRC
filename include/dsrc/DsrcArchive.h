/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCARCHIVE
#define H_DSRCARCHIVE

#include "Globals.h"
#include "Configurable.h"
#include "FastqRecord.h"

namespace dsrc
{

namespace wrap
{

struct FastqRecord;

class DsrcArchive : public Configurable
{
public:
	DsrcArchive();
	~DsrcArchive();

	void StartCompress(const std::string& filename_);
	void WriteNextRecord(const FastqRecord& rec_);
	void FinishCompress();

	void StartDecompress(const std::string& filename_);
	bool ReadNextRecord(FastqRecord& rec_);
	void FinishDecompress();

private:
	struct ArchiveImpl;

	enum DsrcState
	{
		StateNone,
		StateCompression,
		StateDecompression
	};

	DsrcState state;
	std::string filename;

	// internals from DSRC
	ArchiveImpl* impl;

	void FlushChunk();
	bool FeedChunk();

	// hide
	using Configurable::SetStdIoUsing;
	using Configurable::IsStdIoUsing;
	using Configurable::SetThreadsNumber;
	using Configurable::GetThreadsNumber;
};


} // namespace wrap

} // namespace dsrc


#endif // H_DSRCARCHIVE
