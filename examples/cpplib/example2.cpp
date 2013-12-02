#include <Dsrc.h>

#include <iostream>

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
	using namespace dsrc::lib;

	// 1. Open FASTQ file
	//
	FastqFile fastqFile;
	try
	{
		fastqFile.Open(inFilename_);
	}
	catch (const DsrcException& e)
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	// 2. Create and configure DSRC archive
	//
	DsrcArchive archive;
	try
	{
		archive.SetLossyCompression(true);
		archive.SetTagFieldFilterMask(FieldMask().AddField(1).AddField(2).GetMask());
		archive.SetDnaCompressionLevel(2);
		archive.SetQualityCompressionLevel(2);
		archive.SetPlusRepetition(false);								// discard repeated TAG information in "+" lines
		archive.SetFastqBufferSizeMB(256);

		archive.StartCompress(outFilename_);
	}
	catch (const DsrcException& e)
	{
		fastqFile.Close();

		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
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

	std::cerr << "Sucess!\nCompressed records: " << recCount << std::endl;
	return 0;
}


int CompressFile_Solid_Lossless(const char* inFilename_, const char* outFilename_)
{
	using namespace dsrc::lib;

	// 1. Open FASTQ file
	//
	FastqFile fastqFile;
	try
	{
		fastqFile.Open(inFilename_);
	}
	catch (const DsrcException& e)
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	// 2. Create and configure DSRC archive
	//
	DsrcArchive archive;
	try
	{
		archive.SetColorSpace(true);

		archive.SetDnaCompressionLevel(2);
		archive.SetQualityCompressionLevel(1);
		archive.SetPlusRepetition(false);						// discard repeated TAG information in "+" lines
		archive.SetFastqBufferSizeMB(64);

		archive.StartCompress(outFilename_);
	}
	catch (const DsrcException& e)
	{
		fastqFile.Close();

		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
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

	std::cerr << "Sucess!\nCompressed records: " << recCount << std::endl;
	return 0;
}


int DecompressFile(const char* inFilename_, const char* outFilename_)
{
	using namespace dsrc::lib;

	// 1. Open DSRC archive
	//
	DsrcArchive archive;
	try
	{
		archive.StartDecompress(inFilename_);
	}
	catch (const DsrcException& e)
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	// 2. Create output FASTQ file
	//
	FastqFile fastqFile;
	try
	{
		fastqFile.Create(outFilename_);
	}
	catch (const DsrcException& e)
	{
		archive.FinishDecompress();

		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
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

	std::cerr << "Success!\nDecompressed records: " << recCount << std::endl;
	return 0;
}
