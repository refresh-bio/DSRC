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

int CompressFile(const char* inFilename_, const char* outFilename_)
{
	// Define compression configuration settings for DSRC module:
	// lossless mode (by default) with 'fast' performance setting
	// (corresponds to '-m0' switch in DSRC binary)
	DsrcCompressionSettings settings;
	settings.dnaCompressionLevel = 0;
	settings.qualityCompressionLevel = 0;
	settings.fastqBufferSizeMb = 8;

	DsrcModule dsrc;
	if (!dsrc.Compress(inFilename_,
					   outFilename_,
					   settings,		// compress using specified settings
					   2))				// and using 2 processing threads
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << dsrc.GetError() << std::endl;
		return -1;
	}

	std::cerr << "Success!" << std::endl;
	return 0;
}

int CompressFile_lossy(const char* inFilename_, const char* outFilename_)
{
	// Define compression configuration settings for DSRC module:
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

	DsrcModule dsrc;
	if (!dsrc.Compress(inFilename_,
					   outFilename_,
					   settings,		// compress using specified settings
					   2))				// and using 2 processing threads
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << dsrc.GetError() << std::endl;
		return -1;
	}

	std::cerr << "Success!" << std::endl;
	return 0;
}

int DecompressFile(const char* inFilename_, const char* outFilename_)
{
	DsrcModule dsrc;
	if (!dsrc.Decompress(inFilename_,
						 outFilename_,
						 2))			// decompress using 2 processing threads
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << dsrc.GetError() << std::endl;
		return -1;
	}

	std::cerr << "Success!" << std::endl;
	return 0;
}
