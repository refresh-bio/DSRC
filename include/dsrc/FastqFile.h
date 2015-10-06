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

namespace ext
{

class FastqFile
{
public:
	static bool AnalyzeFastqDatasetType(FastqDatasetType& type_, byte* buffer_, uint64 bufferSize_);

	FastqFile();
	~FastqFile();

	bool Open(const std::string& filename_);
	bool Create(const std::string& filename_);
	void Close();

	bool ReadNextRecord(FastqRecord& rec_);
	void WriteNextRecord(const FastqRecord& rec_);

	bool GetFastqDatasetType(FastqDatasetType& type_);

private:
	struct FastqFileImpl;
	FastqFileImpl* fastqImpl;
};


} // namespace ext

} // namespace dsrc

#endif // H_FASTQFILE
