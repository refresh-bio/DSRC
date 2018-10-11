/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "../include/dsrc/Globals.h"
#include "../include/dsrc/DsrcModule.h"

#include <iostream>
#include <cstring>

#include "Common.h"
#include "utils.h"

using namespace dsrc;
using namespace dsrc::ext;
using namespace dsrc::core;

struct InputArguments
{
	enum ModeEnum
	{
		None,
		CompressMode,
		DecompressMode
	};

	static const int MinArguments = 3;

	ModeEnum mode;
	std::string inputFilename;
	std::string outputFilename;
	DsrcCompressionSettings compSettings;

	uint32 threadsNum;
	bool useFastqStdIo;
	uint32 qualityOffset;

	bool verboseMode;

	InputArguments()
		:	mode(None)
		,	threadsNum(1)			// TODO: default
		,	useFastqStdIo(false)
		,	qualityOffset(0)		// TODO: default
		,	verboseMode(false)
	{}
};

void message();
bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_);

int main(int argc_, const char* argv_[])
{
	InputArguments args;
	if (argc_ < InputArguments::MinArguments + 1)
	{
		message();
		return -1;
	}

	if (!parse_arguments(argc_, argv_, args))
	{
		return -1;
	}

	bool success = true;
	DsrcModule dsrc;
	if (args.mode == InputArguments::CompressMode)
	{
		success = dsrc.Compress(args.inputFilename,
								args.outputFilename,
								args.compSettings,
								args.threadsNum,
								args.useFastqStdIo,
								args.qualityOffset);
	}
	else
	{
		success = dsrc.Decompress(args.inputFilename,
								  args.outputFilename,
								  args.threadsNum,
								  args.useFastqStdIo);
	}

	if (!success)
	{
		ASSERT(dsrc.IsError());

		std::cerr << dsrc.GetError();
		return -1;
	}

	if (args.verboseMode && dsrc.GetLog().length() > 0)
	{
		std::cerr << dsrc.GetLog();
	}

	return 0;
}


void message()
{
	std::cerr << "DSRC - DNA Sequence Reads Compressor\n";
	std::cerr << "version: " << DsrcModule::Version() << "\n\n";
	std::cerr << "usage: dsrc <c|d> [options] <input filename> <output filename>\n";
	std::cerr << "compression options:\n";
	std::cerr << "\t-d<n>\t: DNA compression mode: 0-3, default: " << DsrcCompressionSettings::DefaultDnaCompressionLevel << '\n';
	std::cerr << "\t-q<n>\t: Quality compression mode: 0-2, default: " << DsrcCompressionSettings::DefaultQualityCompressionLevel << '\n';
	std::cerr << "\t-f<1,..>: keep only those fields no. in tag field string, default: keep all" << '\n';
	std::cerr << "\t-b<n>\t: FASTQ input buffer size in MB, default: " << DsrcCompressionSettings::DefaultFastqBufferSizeMB << '\n';
	std::cerr << "\t-o<n>\t: Quality offset, default: " << FastqDatasetType::AutoQualityOffsetSelect << " (auto selection)\n";
	std::cerr << "\t-l\t: use Quality lossy mode (Illumina binning scheme), default: false\n";
	std::cerr << "\t-c\t: calculate and check CRC32 checksum calculation per block, default: false\n";

	std::cerr << "automated compression modes:\n";
	std::cerr << "\t-m<n>\t: compression mode, where n:\n";
	std::cerr << "\t * 0\t- fast version with decent ratio (-d0 -q0 -b16)\n";
	std::cerr << "\t * 1\t- slower version with better ratio (-d2 -q2 -b64)\n";
	std::cerr << "\t * 2\t- slow version with best ratio (-d3 -q2 -b256)\n";
	//std::cerr << "\t * 3\t- option (2) with lossy quality and field filtering (-d3 -q2 -b256 -l -f1,2)\n";

	std::cerr << "both compression and decompression options:\n";
	std::cerr << "\t-t<n>\t: processing threads number, default: " << DsrcModule::AvailableHardwareThreadsNum << "(available h/w threads), max: 64" << '\n';
	std::cerr << "\t-s\t: read FASTQ data from stdin (compression) or write FASTQ data to stdout (decompression)\n";
	std::cerr << "\t-v\t: verbose mode, default: false\n\n";

	std::cerr << "usage examples:\n";
	std::cerr << "* compress SRR001471.fastq file saving DSRC archive to SRR001471.dsrc:\n";
	std::cerr << "\tdsrc c SRR001471.fastq SRR001471.dsrc\n";
	std::cerr << "* compress file in the fast mode with CRC32 checking and using 4 threads:\n";
	std::cerr << "\tdsrc c -m0 -c -t4 SRR001471.fastq SRR001471.dsrc\n";
	std::cerr << "* compress file using DNA and Quality compression level 2 and using 512 MB buffer:\n";
	std::cerr << "\tdsrc c -d2 -q2 -b512 SRR001471.fastq SRR001471.dsrc\n";
	std::cerr << "* compress file in the best mode with lossy Quality mode and preserving only 1–4 fields from record IDs:\n";
	std::cerr << "\tdsrc c -m2 -l -f1,2,3,4 SRR001471.fastq SRR001471.dsrc\n";
	std::cerr << "* compress in the best mode reading raw FASTQ data from stdin:\n";
	std::cerr << "\tcat SRR001471.fastq | dsrc c -s -m2 -s SRR001471.dsrc\n";
	std::cerr << "* decompress SRR001471.dsrc archive saving output FASTQ file to SRR001471.out.fastq:\n";
	std::cerr << "\tdsrc d SRR001471.dsrc SRR001471.out.fastq\n";
	std::cerr << "* decompress archive using 4 threads and streaming raw FASTQ data to stdout:\n";
	std::cerr << "\tdsrc d -t4 -s SRR001471.dsrc > SRR001471.out.fastq\n";
}

bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_)
{
	if (argc_ < InputArguments::MinArguments + 1)
		return false;

	if (argv_[1][0] != 'c' && argv_[1][0] != 'd')
	{
		std::cerr << "Error: invalid mode specified\n";
		return false;
	}

	outArgs_.mode = (argv_[1][0] == 'c') ? InputArguments::CompressMode : InputArguments::DecompressMode;

	outArgs_.compSettings = DsrcCompressionSettings::Default();
	DsrcCompressionSettings& compSettings = outArgs_.compSettings;

	outArgs_.threadsNum = DsrcModule::AvailableHardwareThreadsNum;
	outArgs_.qualityOffset = FastqDatasetType::AutoQualityOffsetSelect;
	outArgs_.useFastqStdIo = false;
	outArgs_.verboseMode = false;

	// parse params
	//
	for (int i = 2; i < argc_ - 1; ++i)
	{
		const char* param = argv_[i];
		if (param[0] != '-')
			continue;

		int pval = -1;
		int len = strlen(param);
		if (len > 2)
			pval = to_num((const uchar*)param + 2, len - 2);

		switch (param[1])
		{
			case 't':	outArgs_.threadsNum = pval;						break;
			case 'o':	outArgs_.qualityOffset = pval;					break;
			case 's':	outArgs_.useFastqStdIo = true;					break;
			case 'd':	compSettings.dnaCompressionLevel = pval;		break;
			case 'q':	compSettings.qualityCompressionLevel = pval;	break;
			case 'b':	compSettings.fastqBufferSizeMb = pval;			break;
			case 'l':	compSettings.lossyQualityCompression = true;	break;
			case 'c':	compSettings.calculateCrc32 = true;				break;

			case 'f':
			{
				//pars.fieldsCutFlags = BIT(pval);
				int beg = 2;
				for (int i = 2; i < len; ++i)
				{
					if (param[i] == ',')
					{
						int f = to_num((const uchar*)param + beg, i - beg);
						compSettings.tagPreserveMask |= BIT(f);
						beg = i + 1;
					}
				}
				int f = to_num((const uchar*)param + beg, len - beg);
				compSettings.tagPreserveMask |= BIT(f);

				break;
			}
			case 'm':
			{
				switch (pval)
				{
					//case 3:
					//	pars.tagPreserveFlags = BIT(1) | BIT(2);
					case 2:
						compSettings.dnaCompressionLevel = 3;
						compSettings.qualityCompressionLevel = 2;
						compSettings.fastqBufferSizeMb = 256;
						break;
					case 1:
						compSettings.dnaCompressionLevel = 2;
						compSettings.qualityCompressionLevel = 2;
						compSettings.fastqBufferSizeMb = 64;
						break;
					case 0:
						compSettings.dnaCompressionLevel = 0;
						compSettings.qualityCompressionLevel = 0;
						compSettings.fastqBufferSizeMb = 8;
						break;
				}

				break;
			}

			case 'v':	outArgs_.verboseMode = true;		break;
		}
	}


	// parse input/output file
	//
	outArgs_.inputFilename.clear();
	outArgs_.outputFilename.clear();
	if (!outArgs_.useFastqStdIo)
	{
		outArgs_.inputFilename = argv_[argc_-2];
		outArgs_.outputFilename = argv_[argc_-1];
	}
	else
	{
		if (outArgs_.mode == InputArguments::CompressMode)
			outArgs_.outputFilename = argv_[argc_-1];
		else
			outArgs_.inputFilename = argv_[argc_-1];
	}


	// check params
	//
	if (outArgs_.inputFilename == outArgs_.outputFilename)
	{
		std::cerr << "Error: input and output filenames are the same\n";
		return false;
	}

	//if (outArgs_.verboseMode)
	{
		std::string* fastqFilename = NULL;
		std::string* dsrcFilename = NULL;
		if (outArgs_.mode == InputArguments::CompressMode)
		{
			if (!outArgs_.useFastqStdIo)
				fastqFilename = &outArgs_.inputFilename;
			dsrcFilename = &outArgs_.outputFilename;
		}
		else
		{
			if (!outArgs_.useFastqStdIo)
				fastqFilename = &outArgs_.outputFilename;
			dsrcFilename = &outArgs_.inputFilename;
		}

		if (fastqFilename != NULL && !ends_with(*fastqFilename, ".fastq")
				&& !ends_with(*fastqFilename, ".fq") )
			std::cerr << "Warning: passing a FASTQ file without '.fastq|.fq' extension\n";

		if (dsrcFilename != NULL && !ends_with(*dsrcFilename, ".dsrc"))
			std::cerr << "Warning: passing a DSRC file without '.dsrc' extension\n";
	}

	if (outArgs_.qualityOffset != FastqDatasetType::AutoQualityOffsetSelect
			&& !(outArgs_.qualityOffset >= 33 && outArgs_.qualityOffset <= 64) )
	{
		std::cerr << "Error: invalid Quality offset mode specified [33, 64]\n";
		return false;
	}

	if (outArgs_.threadsNum == 0 || outArgs_.threadsNum > 64)								// limit to 64 ATM ;)
	{
		std::cerr << "Error: invalid thread number specified [1-64]\n";
		return false;
	}

	if (compSettings.dnaCompressionLevel > DsrcCompressionSettings::MaxDnaCompressionLevel)
	{
		std::cerr << "Error: invalid DNA compression mode specified [0-3]\n";
		return false;
	}

	if (compSettings.qualityCompressionLevel > DsrcCompressionSettings::MaxQualityCompressionLevel)
	{
		std::cerr << "Error: invalid Quality compression mode specified [0-2]\n";
		return false;
	}

	if ( !(compSettings.fastqBufferSizeMb >= DsrcCompressionSettings::MinFastqBufferSizeMB
		   && compSettings.fastqBufferSizeMb <= DsrcCompressionSettings::MaxFastqBufferSizeMB) )
	{
		std::cerr << "Error: invalid fastq buffer size specified [1-1024] \n";
		return false;
	}

	return true;
}
