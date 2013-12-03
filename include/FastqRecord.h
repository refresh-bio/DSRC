/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQRECORD
#define H_FASTQRECORD

#include "Globals.h"

namespace dsrc
{

namespace wrap
{

struct FastqRecord
{
	std::string	tag;
	std::string sequence;
	std::string plus;
	std::string quality;
};

} // namespace wrap

} // namespace dsrc


#endif // H_FASTQRECORD
