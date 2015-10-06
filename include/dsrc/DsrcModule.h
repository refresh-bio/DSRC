/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRCMODULE
#define H_DSRCMODULE

#include "Globals.h"

namespace dsrc
{

namespace ext
{


/********
 *
 * DSRC module
 *
 * exposes DSRC compression functionality
 *
 */
class DsrcModule
{
public:
	static const uint32 AvailableHardwareThreadsNum;

	static std::string Version();

	DsrcModule();
	~DsrcModule();

	bool Compress(const std::string& inFastqFilename_,
				  const std::string& outDsrcFilename_,
				  const DsrcCompressionSettings& compSettings_,
				  uint32 threadsNum_,
				  bool useFastqStdIo_ = false,
				  uint32 qualityOffset_ = 0);

	bool Decompress(const std::string& inDsrcFilename_,
					const std::string& outFastqFilename_,
					uint32 threadsNum_,
					bool useFastqStdIo_ = false);

	bool IsError() const;
	const std::string& GetError() const;
	void ClearError();

	const std::string& GetLog() const;
	void ClearLog();

private:
	std::string errorMsg;
	std::string logMsg;

	void AddError(const std::string& err_);
	void SetError(const std::string& err_);
	void AddLog(const std::string& log_);
};

} // namespace ext

} // namespace dsrc


#endif // H_DSRCMODULE
