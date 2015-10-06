#include <iostream>
#include <Dsrc.h>

using namespace dsrc::lib;


void PrintMessage();
int CompressFile_Illumina_Lossy(const char* inFilename_, const char* outFilename_);
int CompressFile_Solid_Lossless(const char* inFilename_, const char* outFilename_);
int DecompressFile(const char* inFilename_, const char* outFilename_);


int main(int argc_, char* argv_[])
{
	if ( !(argc_ == 4 || argc_ == 5) || !(argv_[1][0] == 'c' || argv_[1][0] == 'd') )
	{
		PrintMessage();
		return -1;
	}

	bool useSolid = argc_ == 5 && argv_[2][0] == 'S';
	char option = argv_[1][0];
	const char* inFilename = argv_[argc_-2];
	const char* outFilename = argv_[argc_-1];

	if (option == 'c')
	{
		if (useSolid)
			return CompressFile_Solid_Lossless(inFilename, outFilename);
		else
			return CompressFile_Illumina_Lossy(inFilename, outFilename);
	}
	if (option == 'd')
	{
		return DecompressFile(inFilename, outFilename);
	}

	return -1;
}


void PrintMessage()
{
	std::cerr << "usage: example1 <c|d> [S] <in_filename> <out_filename>" << std::endl;
}


int CompressFile_Illumina_Lossy(const char* inFilename_, const char* outFilename_)
{
	// 1. Open FASTQ file
	//
	FastqFile fastqFile;
	if (!fastqFile.Open(inFilename_))		// TODO: analyze the dataset type
	{
		std::cerr << "Error: cannot open FASTQ file: " << inFilename_ << std::endl;
		return -1;
	}


	// 2. Create and configure DSRC archive records writer
	//
	// specify DSRC archive compression settings
	DsrcCompressionSettings settings;
	settings.lossyQualityCompression = true;
	settings.tagPreserveMask = FieldMask().AddField(1).AddField(2).GetMask();
	settings.dnaCompressionLevel = 2;
	settings.qualityCompressionLevel = 2;
	settings.fastqBufferSizeMb = 256;

	// specify explicitely input FASTQ dataset type
	FastqDatasetType datasetType;
	datasetType.plusRepetition = false;
	datasetType.colorSpace = false;

	// initialize compression routine with prepared settings
	DsrcArchiveRecordsWriter archive;
	if (!archive.StartCompress(outFilename_,
							   settings,
							   datasetType))
	{
		fastqFile.Close();

		std::cerr << "Error!" << std::endl;
		std::cerr << archive.GetError() << std::endl;
		return -1;
	}


	// 3. When configured, start compression
	//
	int recCount = 0;
	FastqRecord rec;
	while (fastqFile.ReadNextRecord(rec))
	{
		archive.WriteNextRecord(rec);
		recCount++;
	}

	// 4. Finish compression and close FASTQ file
	//
	archive.FinishCompress();

	fastqFile.Close();

	std::cerr << "Sucess!" << std::endl;
	std::cerr << "Compressed records: " << recCount << std::endl;
	return 0;
}


int CompressFile_Solid_Lossless(const char* inFilename_, const char* outFilename_)
{
	using namespace dsrc::lib;

	// 1. Open FASTQ file
	//
	FastqFile fastqFile;
	if (!fastqFile.Open(inFilename_))			// TODO: analyze the dataset type
	{
		std::cerr << "Error: cannot open FASTQ file: " << inFilename_ << std::endl;
		return -1;
	}

	// 2. Create and configure DSRC archive
	//
	DsrcCompressionSettings settings;
	settings.dnaCompressionLevel = 2;
	settings.qualityCompressionLevel = 1;
	settings.fastqBufferSizeMb = 64;

	FastqDatasetType datasetType;
	datasetType.colorSpace = true;
	datasetType.plusRepetition = false;	// discard repeated TAG information in "+" lines

	DsrcArchiveRecordsWriter archive;
	if (!archive.StartCompress(outFilename_,
							   settings,
							   datasetType))
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << archive.GetError() << std::endl;
		return -1;
	}

	// 3. When configured, start compression
	//
	long int recCount = 0;
	FastqRecord rec;
	while (fastqFile.ReadNextRecord(rec))
	{
		archive.WriteNextRecord(rec);
		recCount++;
	}

	// 4. Finish compression and close FASTQ file
	//
	archive.FinishCompress();

	fastqFile.Close();

	std::cerr << "Sucess!" << std::endl;
	std::cerr << "Compressed records: " << recCount << std::endl;
	return 0;
}


int DecompressFile(const char* inFilename_, const char* outFilename_)
{
	// 1. Open DSRC archive
	//
	DsrcArchiveRecordsReader archive;
	if (!archive.StartDecompress(inFilename_))
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << archive.GetError() << std::endl;
		return -1;
	}

	// 2. Create output FASTQ file
	//
	FastqFile fastqFile;
	if (!fastqFile.Create(outFilename_))
	{
		archive.FinishDecompress();

		std::cerr << "Error: cannot open FASTQ file: " << inFilename_ << std::endl;
		return -1;
	}

	// 3. Start decompression
	//
	long int recCount = 0;
	FastqRecord rec;
	while (archive.ReadNextRecord(rec))
	{
		fastqFile.WriteNextRecord(rec);
		recCount++;
	}

	// 4. Finish decompression and close FASTQ file
	//
	fastqFile.Close();

	archive.FinishDecompress();

	std::cerr << "Success!" << std::endl;
	std::cerr << "Decompressed records: " << recCount << std::endl;
	return 0;
}
