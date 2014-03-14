/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "../include/dsrc/Globals.h"

#include <iostream>
#include <cstring>

#include "DsrcOperator.h"
#include "utils.h"

using namespace dsrc;
using namespace dsrc::core;
using namespace dsrc::comp;

struct InputArguments
{
	enum ModeEnum
	{
		CompressMode,
		DecompressMode
	};

	static const int MinArguments = 3;

	ModeEnum mode;
	InputParameters params;
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

	IDsrcOperator* op = NULL;
	if (args.params.threadNum == 1)
	{
		if (args.mode == InputArguments::CompressMode)
			op = new DsrcCompressorST();
		else
			op = new DsrcDecompressorST();
	}
	else
	{
		if (args.mode == InputArguments::CompressMode)
			op = new DsrcCompressorMT();
		else
			op = new DsrcDecompressorMT();
	}

	if (!op->Process(args.params))
	{
		ASSERT(op->IsError());
		std::cerr << op->GetError();
		delete op;
		return -1;
	}

	delete op;
	return 0;
}


void message()
{
	std::cerr << "DSRC - DNA Sequence Reads Compressor\n";
	std::cerr << "version: 2.00 RC @ 01.12.2013\n\n";
	std::cerr << "usage: dsrc <c|d> [options] <input filename> <output filename>\n";
	std::cerr << "compression options:\n";
	std::cerr << "\t-d<n>\t: DNA compression mode: 0-3, default: " << InputParameters::DefaultDnaCompressionLevel << '\n';
	std::cerr << "\t-q<n>\t: Quality compression mode: 0-2, default: " << InputParameters::DefaultQualityCompressionLevel << '\n';
	std::cerr << "\t-f<1,..>: keep only those fields no. in tag field string, default: keep all" << '\n';
	std::cerr << "\t-b<n>\t: FASTQ input buffer size in MB, default: " << InputParameters::DefaultFastqBufferSizeMB << '\n';
	std::cerr << "\t-o<n>\t: Quality offset, default: " << InputParameters::DefaultQualityOffset << '\n';
	std::cerr << "\t-l\t: use Quality lossy mode (Illumina binning scheme), default: " << InputParameters::DefaultLossyCompressionMode << '\n';
	std::cerr << "\t-c\t: calculate and check CRC32 checksum calculation per block, default: " << InputParameters::DefaultCalculateCrc32 << '\n';

	std::cerr << "automated compression modes:\n";
	std::cerr << "\t-m<n>\t: compression mode, where n:\n";
	std::cerr << "\t * 0\t- fast version with decent ratio (-d0 -q0 -b16)\n";
	std::cerr << "\t * 1\t- slower version with better ratio (-d2 -q2 -b64)\n";
	std::cerr << "\t * 2\t- slow version with best ratio (-d3 -q2 -b256)\n";
	//std::cerr << "\t * 3\t- option (2) with lossy quality and field filtering (-d3 -q2 -b256 -l -f1,2)\n";

	std::cerr << "both compression and decompression options:\n";
	std::cerr << "\t-t<n>\t: processing threads number, default (available h/w threads): " << IDsrcOperator::AvailableHardwareThreadsNum << ", max: 64" << '\n';
	std::cerr << "\t-s\t: use stdin/stdout for reading/writing FASTQ data\n";
}

bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_)
{
	if (argc_ < InputArguments::MinArguments + 1)
		return false;

	if (argv_[1][0] != 'c' && argv_[1][0] != 'd')
	{
		std::cerr << "Error: invalid Quality offset mode specified\n";
	}

	outArgs_.mode = (argv_[1][0] == 'c') ? InputArguments::CompressMode : InputArguments::DecompressMode;

	outArgs_.params = InputParameters::Default();
	InputParameters& pars = outArgs_.params;
	pars.threadNum = IDsrcOperator::AvailableHardwareThreadsNum;

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
			case 'o':	pars.qualityOffset = pval;			break;
			case 'd':	pars.dnaCompressionLevel = pval;	break;
			case 'q':	pars.qualityCompressionLevel = pval; break;
			case 't':	pars.threadNum = pval;				break;
			case 'b':	pars.fastqBufferSizeMB = pval;		break;
			case 'l':	pars.lossyCompression = true;		break;
			case 'c':	pars.calculateCrc32 = true;			break;
			case 's':	pars.useFastqStdIo = true;			break;
			case 'f':
			{
				//pars.fieldsCutFlags = BIT(pval);
				int beg = 2;
				for (int i = 2; i < len; ++i)
				{
					if (param[i] == ',')
					{
						int f = to_num((const uchar*)param + beg, i - beg);
						pars.tagPreserveFlags |= BIT(f);
						beg = i + 1;
					}
				}
				int f = to_num((const uchar*)param + beg, len - beg);
				pars.tagPreserveFlags |= BIT(f);

				break;
			}
			case 'm':
			{
				switch (pval)
				{
					//case 3:
					//	pars.tagPreserveFlags = BIT(1) | BIT(2);
					case 2:
						pars.dnaCompressionLevel = 3;
						pars.qualityCompressionLevel = 2;
						pars.fastqBufferSizeMB = 256;
						break;
					case 1:
						pars.dnaCompressionLevel = 2;
						pars.qualityCompressionLevel = 1;
						pars.fastqBufferSizeMB = 64;
						break;
					case 0:
						pars.dnaCompressionLevel = 0;
						pars.qualityCompressionLevel = 0;
						pars.fastqBufferSizeMB = 8;
						break;
				}

				break;
			}
		}
	}


	// check params
	//
	if (!pars.useFastqStdIo)
	{
		pars.inputFilename = argv_[argc_-2];
		pars.outputFilename = argv_[argc_-1];
	}
	else
	{
		if (outArgs_.mode == InputArguments::CompressMode)
			pars.outputFilename = argv_[argc_-1];
		else
			pars.inputFilename = argv_[argc_-1];
	}

	if (pars.qualityOffset != fq::FastqDatasetType::AutoQualityOffset && !(pars.qualityOffset >= 33 && pars.qualityOffset <= 64) )
	{
		std::cerr << "Error: invalid Quality offset mode specified [33, 64]\n";
		return false;
	}

	if (pars.dnaCompressionLevel > 3)
	{
		std::cerr << "Error: invalid DNA compression mode specified [0-3]\n";
		return false;
	}

	if (pars.qualityCompressionLevel > 2)
	{
		std::cerr << "Error: invalid Quality compression mode specified [0-2]\n";
		return false;
	}

	if (pars.threadNum == 0 || pars.threadNum > 64)
	{
		std::cerr << "Error: invalid thread number specified [1-64]\n";
		return false;
	}

	if ( !(pars.fastqBufferSizeMB >= 1 && pars.fastqBufferSizeMB <= 1024) )		// limit to 1GB ATM ;)
	//	   && (pars.fastqBufferSizeMB & (pars.fastqBufferSizeMB - 1)) == 0) )
	{
		std::cerr << "Error: invalid fastq buffer size specified [1-1024] \n";
		return false;
	}

	return true;
}
