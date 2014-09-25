/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "../include/dsrc/Globals.h"

#include <iostream>

#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
namespace th = boost;
#else
#include <thread>
namespace th = std;
#endif

#include "../include/dsrc/DsrcModule.h"

#include "DsrcOperator.h"


namespace dsrc
{

using namespace comp;

namespace wrap
{

const uint32 DsrcModule::HardwareThreadsNo = th::thread::hardware_concurrency();

DsrcModule::DsrcModule()
{
	SetThreadsNumber(HardwareThreadsNo);
}

DsrcModule::~DsrcModule()
{

}

void DsrcModule::Compress(const std::string &inputFilename_, const std::string &outputFilename_)
{
	IDsrcOperator* dsrc;

	if (GetThreadsNumber() > 0)
		dsrc = new DsrcCompressorMT();
	else
		dsrc = new DsrcCompressorST();

	InputParameters params = *(const InputParameters*)GetInputParameters();
	params.inputFilename = inputFilename_;
	params.outputFilename = outputFilename_;
	if (!dsrc->Process(params))
	{
		std::string err = dsrc->GetError();
		delete dsrc;

		throw DsrcException(err);
	}

	delete dsrc;
}

void DsrcModule::Decompress(const std::string &inputFilename_, const std::string &outputFilename_)
{
	IDsrcOperator* dsrc;

	if (GetThreadsNumber() > 0)
		dsrc = new DsrcDecompressorMT();
	else
		dsrc = new DsrcDecompressorST();

	InputParameters params = *(const InputParameters*)GetInputParameters();
	params.inputFilename = inputFilename_;
	params.outputFilename = outputFilename_;
	if (!dsrc->Process(params))
	{
		std::string err = dsrc->GetError();
		delete dsrc;

		throw DsrcException(err);
	}

	delete dsrc;
}

void DsrcModule::Usage()
{
	std::cerr << "DSRC 2.00b release\n";
}

} // namespace wrap

} // namespace dsrc
