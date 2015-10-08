/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCARCHIVE
#define H_DSRCARCHIVE

#include "Globals.h"
#include "FastqRecord.h"

namespace dsrc
{

namespace ext
{

struct FastqRecord;

/********
 *
 * DSRC archive writer/reader interfaces
 *
 */
class DsrcArchiveWriter
{
public:
	DsrcArchiveWriter();
	virtual ~DsrcArchiveWriter();

	virtual bool StartCompress(const std::string& dsrcFilename_,
							   const DsrcCompressionSettings& compressionSettings_,
							   uint32 threadsNum_ = 1,				// for future use, currently only 1
							   uint32 qualityOffset_ = 0) = 0;		// auto selection of quality offset: 0
	virtual void FinishCompress() = 0;

	bool IsError() const;
	const std::string& GetError() const;
	void ClearError();

protected:
	struct ArchiveWriterImpl;
	ArchiveWriterImpl* archiveImpl;
};


class DsrcArchiveReader
{
public:
	DsrcArchiveReader();
	virtual ~DsrcArchiveReader();

	virtual bool StartDecompress(const std::string& dsrcFilename_,
								 uint32 threadsNum_ = 1) = 0;		// for future use, currently only 1
	virtual void FinishDecompress() = 0;

	bool GetCompressionSettings(DsrcCompressionSettings& settings_) const;

	bool IsError() const;
	const std::string& GetError() const;
	void ClearError();

protected:
	struct ArchiveReaderImpl;
	ArchiveReaderImpl* archiveImpl;
};


/********
 *
 * DSRC archive records writer/reader classes
 *
 * implements writing to / reading from DSRC archives
 * record by record
 *
 */
class DsrcArchiveRecordsWriter : public DsrcArchiveWriter
{
public:
	DsrcArchiveRecordsWriter();
	~DsrcArchiveRecordsWriter();

	bool StartCompress(const std::string& dsrcFilename_,
					   const DsrcCompressionSettings& compressionSettings_,
					   uint32 threadsNum_ = 1,				// for future use, currently only 1
					   uint32 qualityOffset_ = 0);			// auto selection of quality offset: 0
	void FinishCompress();

	bool WriteNextRecord(const FastqRecord& rec_);

protected:
	struct RecordsWriterImpl;
	RecordsWriterImpl* writerImpl;
};


class DsrcArchiveRecordsReader : public DsrcArchiveReader
{
public:
	DsrcArchiveRecordsReader();
	~DsrcArchiveRecordsReader();

	bool StartDecompress(const std::string &dsrcFilename_,
						 uint32 threadsNum_ = 1);
	void FinishDecompress();

	bool ReadNextRecord(FastqRecord& rec_);

protected:
	struct RecordsReaderImpl;
	RecordsReaderImpl* readerImpl;
};


/********
 *
 * DSRC archive block writer/reader classes
 *
 * implements a single threaded (at the moment)
 * DSRC archive writer/reader block by block
 *
 */
class DsrcArchiveBlocksWriterST : public DsrcArchiveWriter
{
public:
	DsrcArchiveBlocksWriterST();
	~DsrcArchiveBlocksWriterST();

	bool StartCompress(const std::string& dsrcFilename_,
					   const DsrcCompressionSettings& compressionSettings_,
					   uint32 threadsNum_ = 1,				// for future use, currently only 1
					   uint32 qualityOffset_ = 0);			// auto selection of quality offset: 0
	void FinishCompress();

	// returns the number of compressed bytes written
	unsigned long WriteNextBlock(const char* buffer_, unsigned long bufferSize_);

private:
	struct BlockWriterImpl;
	BlockWriterImpl* writerImpl;
};


class DsrcArchiveBlocksReaderST : public DsrcArchiveReader
{
public:
	DsrcArchiveBlocksReaderST();
	~DsrcArchiveBlocksReaderST();

	bool StartDecompress(const std::string &dsrcFilename_,
						 uint32 threadsNum_ = 1);			// for future use, currently only 1
	void FinishDecompress();

	// returns number of bytes read
	unsigned long ReadNextBlock(char* buffer_, unsigned long bufferSize_);

private:
	struct BlockReaderImpl;
	BlockReaderImpl* readerImpl;
};


typedef DsrcArchiveBlocksWriterST DsrcArchiveBlocksWriter;
typedef DsrcArchiveBlocksReaderST DsrcArchiveBlocksReader;


} // namespace ext

} // namespace dsrc


#endif // H_DSRCARCHIVE
