/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "StdStream.h"

#include <stdio.h>

namespace dsrc
{

namespace core
{


// ********************************************************************************************
int64 StdStreamReader::Read(uchar *mem_, uint64 size_)
{
	int64 n = fread(mem_, 1, size_, stdin);
	return n;
}

// ********************************************************************************************
int64 StdStreamWriter::Write(const uchar *mem_, uint64 size_)
{
	int64 n = fwrite(mem_, 1, size_, stdout);
	return n;
}


} // namespace core

} // namespace dsrc
