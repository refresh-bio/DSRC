#include <Dsrc.h>

#include <iostream>

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
	using namespace dsrc::lib;


	// Configure DSRC compressor
	//
	try
	{
		DsrcModule dsrc;
		dsrc.SetLossyCompression(true);
		dsrc.SetTagFieldFilterMask(FieldMask().AddField(1).AddField(2).GetMask());
		dsrc.SetDnaCompressionLevel(2);
		dsrc.SetQualityCompressionLevel(2);

		dsrc.SetFastqBufferSizeMB(256);
		dsrc.SetThreadsNumber(2);

		// additionally perform CRC32 checking when compressing
		dsrc.SetCrc32Checking(true);

		dsrc.Compress(inFilename_, outFilename_);
	}
	catch (const DsrcException& e)
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	std::cerr << "Success!" << std::endl;
	return 0;
}


int DecompressFile(const char* inFilename_, const char* outFilename_)
{
	using namespace dsrc::lib;

	try
	{
		DsrcModule dsrc;
		dsrc.Decompress(inFilename_, outFilename_);
	}
	catch (const DsrcException& e)
	{
		std::cerr << "Error!" << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	std::cerr << "Success!" << std::endl;
	return 0;
}
