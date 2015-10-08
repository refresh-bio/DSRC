#include <iostream>
#include <Dsrc.h>

using namespace dsrc::lib;


void PrintMessage();
int CompressFile(const char* inFilename_, const char* outFilename_);
int CompressFile_lossy(const char* inFilename_, const char* outFilename_);
int DecompressFile(const char* inFilename_, const char* outFilename_);


int main(int argc_, char* argv_[])
{
	if ( !(argc_ == 4 || argc_ == 5) || !(argv_[1][0] == 'c' || argv_[1][0] == 'd') )
	{
		PrintMessage();
		return -1;
	}

	char option = argv_[1][0];
	bool lossyMode = argc_ == 5 && argv_[2][0] == 'L';
	const char* inFilename = argv_[argc_-2];
	const char* outFilename = argv_[argc_-1];

	if (option == 'c')
	{
		if (lossyMode)
			return CompressFile_lossy(inFilename, outFilename);
		else
			return CompressFile(inFilename, outFilename);
	}
	else if (option == 'd')
	{
		return DecompressFile(inFilename, outFilename);
	}

	return -1;
}

void PrintMessage()
{
	std::cerr << "usage: example1 <c|d> [L] <in_filename> <out_filename>" << std::endl;
}


// RAII buffer wrapper (but better use smart pointers, if possible)
struct FastqBuffer
{
	FastqBuffer(unsigned long int capacity_)
		:	capacity(capacity_)
	{
		data = new char[capacity_];
	}

	~FastqBuffer()
	{
		delete data;
	}

	const unsigned long int capacity;
	char* data;
};


int CompressFile(const char* inFilename_, const char* outFilename_)
{
	// 1. Open FASTQ file
	//
	FastqFileBlocksReader fastqFile;
	if (!fastqFile.Open(inFilename_))		// TODO: analyze the dataset type
	{
		std::cerr << "Error: cannot open FASTQ file: " << inFilename_ << std::endl;
		return -1;
	}

	// 2. Create and configure DSRC archive records writer with comp. settings:
	// lossless mode (by default) with 'fast' performance setting
	// (corresponds to '-m0' switch in DSRC binary)
	DsrcCompressionSettings settings;
	settings.dnaCompressionLevel = 0;
	settings.qualityCompressionLevel = 0;
	settings.fastqBufferSizeMb = 8;

	// initialize compression routine with prepared settings
	DsrcArchiveBlocksWriter archive;
	if (!archive.StartCompress(outFilename_, settings))
	{
		fastqFile.Close();

		std::cerr << "Error!" << std::endl;
		std::cerr << archive.GetError() << std::endl;
		return -1;
	}

	// 3. When configured, start compression
	//
	FastqBuffer buffer(settings.fastqBufferSizeMb << 20);
	int blocksCount = 0;
	do
	{
		unsigned long int blockSize = fastqFile.ReadNextBlock(buffer.data, buffer.capacity);

		if (blockSize == 0)
			break;

		archive.WriteNextBlock(buffer.data, blockSize);
		blocksCount++;
	}
	while (1);

	// 4. Finish compression and close FASTQ file
	//
	archive.FinishCompress();

	fastqFile.Close();

	std::cerr << "Sucess!" << std::endl;
	std::cerr << "Compressed blocks: " << blocksCount << std::endl;
	return 0;
}


int CompressFile_lossy(const char* inFilename_, const char* outFilename_)
{
	// 1. Open FASTQ file
	//
	FastqFileBlocksReader fastqFile;
	if (!fastqFile.Open(inFilename_))		// TODO: analyze the dataset type
	{
		std::cerr << "Error: cannot open FASTQ file: " << inFilename_ << std::endl;
		return -1;
	}

	// 2. Create and configure DSRC archive records writer with settings:
	// use 'max' ratio compression mode
	// (corresponds to '-m2' switch in DSRC binary)
	DsrcCompressionSettings settings;
	settings.dnaCompressionLevel = 3;
	settings.qualityCompressionLevel = 2;
	settings.fastqBufferSizeMb = 256;

	// set lossy mode with Illumina 8-binning qualities scheme
	settings.lossyQualityCompression = true;

	// keep only 2 first tokens in read tag field
	settings.tagPreserveMask = FieldMask().AddField(1).AddField(2).GetMask();

	// initialize compression routine with prepared settings
	DsrcArchiveBlocksWriter archive;
	if (!archive.StartCompress(outFilename_, settings))
	{
		fastqFile.Close();

		std::cerr << "Error!" << std::endl;
		std::cerr << archive.GetError() << std::endl;
		return -1;
	}

	// 3. When configured, start compression
	//
	FastqBuffer buffer(settings.fastqBufferSizeMb << 20);
	int blocksCount = 0;
	do
	{
		unsigned long int blockSize = fastqFile.ReadNextBlock(buffer.data, buffer.capacity);

		if (blockSize == 0)
			break;

		archive.WriteNextBlock(buffer.data, blockSize);
		blocksCount++;
	}
	while (1);

	// 4. Finish compression and close FASTQ file
	//
	archive.FinishCompress();

	fastqFile.Close();

	std::cerr << "Sucess!" << std::endl;
	std::cerr << "Compressed blocks: " << blocksCount << std::endl;
	return 0;
}


int DecompressFile(const char* inFilename_, const char* outFilename_)
{
	// 1. Open DSRC archive
	//
	DsrcArchiveBlocksReader archive;
	if (!archive.StartDecompress(inFilename_))
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << archive.GetError() << std::endl;
		return -1;
	}

	// 2. Create output FASTQ file
	//
	FastqFileBlocksWriter fastqFile;
	if (!fastqFile.Open(outFilename_))
	{
		archive.FinishDecompress();

		std::cerr << "Error: cannot open FASTQ file: " << inFilename_ << std::endl;
		return -1;
	}

	// 3. Start decompression
	//
	// get the block size of the initial FASTQ file
	// for better performance and memory utilisation
	DsrcCompressionSettings settings;
	archive.GetCompressionSettings(settings);
	FastqBuffer buffer(settings.fastqBufferSizeMb << 20);

	long int blocksCount = 0;
	do
	{
		unsigned long blockSize = archive.ReadNextBlock(buffer.data, buffer.capacity);

		if (blockSize == 0)
			break;

		fastqFile.WriteNextBlock(buffer.data, blockSize);
		blocksCount++;
	}
	while (1);

	// 4. Finish decompression and close FASTQ file
	//
	fastqFile.Close();

	archive.FinishDecompress();

	std::cerr << "Success!" << std::endl;
	std::cerr << "Decompressed blocks: " << blocksCount << std::endl;
	return 0;
}
