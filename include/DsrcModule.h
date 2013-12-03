/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCMODULE
#define H_DSRCMODULE

#include "Globals.h"
#include "Configurable.h"

namespace dsrc
{

namespace wrap
{

class DsrcModule : public Configurable
{
public:
	DsrcModule();
	~DsrcModule();

	void Compress(const std::string& inputFilename_, const std::string& outputFilename_);
	void Decompress(const std::string& inputFilename_, const std::string& outputFilename_);
	void Usage();

private:
	static const uint32 HardwareThreadsNo;

	// hide
	using Configurable::IsColorSpace;
	using Configurable::SetColorSpace;
	using Configurable::IsPlusRepetition;
	using Configurable::SetPlusRepetition;
};

} // namespace wrap

} // namespace dsrc


#endif // H_DSRCMODULE
