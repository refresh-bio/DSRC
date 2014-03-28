/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_ERRORHANDLER
#define H_ERRORHANDLER

#include "../include/dsrc/Globals.h"

#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
namespace th = boost;
#else
#include <mutex>
#include <atomic>
namespace th = std;
#endif


namespace dsrc
{

namespace core
{

class ErrorHandler
{
public:
	virtual ~ErrorHandler() {}

	virtual bool IsError()
	{
		return isError;
	}

	virtual std::string GetError()
	{
		return error;
	}

	virtual void SetError(const std::string& err_)
	{
		error = err_;
		isError = true;
	}

protected:
	bool isError;
	std::string error;
};

class MultithreadedErrorHandler : public ErrorHandler
{
public:
	bool IsError()
	{
		th::lock_guard<th::mutex> lock(mutex);
		return isError;
	}

	std::string GetError()
	{
		th::lock_guard<th::mutex> lock(mutex);
		return error;
	}

	void SetError(const std::string& err_)
	{
		th::lock_guard<th::mutex> lock(mutex);
		if (isError)
			return;
		error = err_;
		isError = true;
	}

private:
	th::mutex mutex;
};

} // namespace core

} // namespace dsrc

#endif // H_ERRORHANDLER
