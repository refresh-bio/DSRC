#include <iostream>
#include <Dsrc.h>

using namespace dsrc::lib;


void PrintMessage();
int CompressFile(const char* inFilename_, const char* outFilename_);
int DecompressFile(const char* inFilename_, const char* outFilename_);

int main(int argc_, char* argv_[])
{
	if ( argc_ != 4 || !(argv_[1][0] == 'c' || argv_[1][0] == 'd') )
	{
		PrintMessage();
		return -1;
	}

	char option = argv_[1][0];
	const char* inFilename = argv_[argc_-2];
	const char* outFilename = argv_[argc_-1];

	if (option == 'c')
		return CompressFile(inFilename, outFilename);

	if (option == 'd')
		return DecompressFile(inFilename, outFilename);

	return -1;
}


void PrintMessage()
{
	std::cerr << "usage: example1 <c|d> <in_filename> <out_filename>" << std::endl;
}


int CompressFile(const char* inFilename_, const char* outFilename_)
{
	// Configure DSRC compressor
	//
	DsrcCompressionSettings settings;
	settings.lossyQualityCompression = true;
	settings.tagPreserveMask = FieldMask().AddField(1).AddField(2).GetMask();
	settings.dnaCompressionLevel = 3;
	settings.dnaCompressionLevel = 2;
	settings.fastqBufferSizeMb = 256;

	// additionally compute and perform CRC32 checking when compressing
	settings.calculateCrc32 = true;

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
