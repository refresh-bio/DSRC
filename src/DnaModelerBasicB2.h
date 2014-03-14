/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DNAMODELERBASICB2
#define H_DNAMODELERBASICB2

#include "../include/dsrc/Globals.h"

#include "DnaModeler.h"
#include "Stats.h"
#include "Fastq.h"
#include "BitMemory.h"

namespace dsrc
{

namespace comp
{

class DnaModelerBasicB2 : public IDnaModeler
{
public:
	void ProcessStats(const DnaStats& stats_)
	{
		ASSERT(stats_.symbolCount <= 4);
	}

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_)
	{
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			const fq::FastqRecord& r = records_[i];
			for (uint32 j = 0; j < r.sequenceLen; ++j)
			{
				ASSERT(r.sequence[j] < 4);
				writer_.Put2Bits(r.sequence[j]);
			}
		}
		writer_.FlushPartialWordBuffer();
	}

	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_)
	{
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			fq::FastqRecord& r = records_[i];
			for (uint32 j = 0; j < r.sequenceLen; ++j)
			{
				r.sequence[j] = reader_.Get2Bits();
				ASSERT(r.sequence[j] < 4);
			}
		}
		reader_.FlushInputWordBuffer();
	}
};

} // namespace comp

} // namespace dsrc

#endif // H_DNAMODELERBASICB2
