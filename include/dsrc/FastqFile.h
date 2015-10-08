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

/********
 *
 * FASTQ file writer/reader interfaces
 *
 */
class FastqFileReader
{
public:
	static bool AnalyzeFastqDatasetType(FastqDatasetType& type_, byte* buffer_, uint64 bufferSize_);

	FastqFileReader();
	virtual ~FastqFileReader();

	virtual bool Open(const std::string& filename_) = 0;
	virtual void Close() = 0;
};


class FastqFileWriter
{
public:
	FastqFileWriter();
	virtual ~FastqFileWriter();

	virtual bool Open(const std::string& filename_) = 0;
	virtual void Close() = 0;
};


/********
 *
 * FASTQ file records writer/reader classes
 *
 * implements writing to / reading from FASTQ files
 * record by record
 *
 */
class FastqFileRecordsReader : public FastqFileReader
{
public:
	FastqFileRecordsReader();
	~FastqFileRecordsReader();

	bool Open(const std::string& filename_);
	void Close();

	bool ReadNextRecord(FastqRecord& rec_);

private:
	struct FastqRecordsReaderImpl;
	FastqRecordsReaderImpl* readerImpl;
};


class FastqFileRecordsWriter : public FastqFileWriter
{
public:
	FastqFileRecordsWriter();
	~FastqFileRecordsWriter();

	bool Open(const std::string& filename_);
	void Close();

	bool WriteNextRecord(const FastqRecord& rec_);

private:
	struct FastqRecordsWriterImpl;
	FastqRecordsWriterImpl* writerImpl;
};


/********
 *
 * FASTQ file block writer/reader classes
 *
 * implements a single threaded
 * FASTQ file writer/reader block by block
 *
 */
class FastqFileBlocksReader : public FastqFileReader
{
public:
	static const uint64 MinBufferSize = 32 * 1024;

	FastqFileBlocksReader();
	~FastqFileBlocksReader();

	bool Open(const std::string& filename_);
	void Close();

	unsigned long int ReadNextBlock(char* buffer_, unsigned long int bufferSize_);

private:
	struct FastqBlocksReaderImpl;
	FastqBlocksReaderImpl* readerImpl;
};


class FastqFileBlocksWriter : public FastqFileWriter
{
public:
	FastqFileBlocksWriter();
	~FastqFileBlocksWriter();

	bool Open(const std::string& filename_);
	void Close();

	unsigned long int WriteNextBlock(char* buffer_, unsigned long int dataSize_);

private:
	struct FastqBlocksWriterImpl;
	FastqBlocksWriterImpl* writerImpl;
};


} // namespace ext

} // namespace dsrc

#endif // H_FASTQFILE
