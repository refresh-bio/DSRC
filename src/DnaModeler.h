/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DNAMODELER
#define H_DNAMODELER

#include "../include/dsrc/Globals.h"

#include "Common.h"

namespace dsrc
{

namespace comp
{

class IDnaModeler
{
public:
	virtual ~IDnaModeler() {}

	virtual void ProcessStats(const DnaStats& stats_) = 0;

	virtual void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_) = 0;
	virtual void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_) = 0;
};

} // namespace comp

} // namespace dsrc


#endif // H_DNAMODELER
