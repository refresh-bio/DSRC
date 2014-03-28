/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQ_SRC
#define H_FASTQ_SRC

#include "../include/dsrc/Globals.h"
#include "Common.h"

#include <vector>

#include "Buffer.h"

#include "utils.h"

namespace dsrc
{

namespace fq
{

typedef core::DataChunk FastqDataChunk;

// 8B padding
struct FastqRecord
{
	uchar *title;
	uchar *sequence;
	uchar *quality;

	uint16 titleLen;
	uint16 sequenceLen;		// can be specialized
	uint16 qualityLen;		// can be specialized
	uint16 truncatedLen;	// can be specialized

	FastqRecord()
		:	title(NULL)
		,	sequence(NULL)
		,	quality(NULL)
		,	titleLen(0)
		,	sequenceLen(0)
		,	qualityLen(0)
		,	truncatedLen(0)
	{}

	void Reset()
	{
		title = NULL;
		titleLen = 0;
		sequence = NULL;
		sequenceLen = 0;
		quality = NULL;
		qualityLen = 0;
		truncatedLen = 0;
	}
};


struct FastqChecksum
{
	enum CheksumFlags
	{
		CALC_TAG		= BIT(0),
		CALC_SEQUENCE	= BIT(1),
		CALC_QUALITY	= BIT(2),
		CALC_NONE		= 0,
		CALC_ALL		= CALC_TAG | CALC_SEQUENCE | CALC_QUALITY
	};

	uint32 tag;
	uint32 sequence;
	uint32 quality;

	FastqChecksum()
		:	tag(0)
		,	sequence(0)
		,	quality(0)
	{}

	void Reset()
	{
		tag = 0;
		sequence = 0;
		quality = 0;
	}
};

} // namespace fq

} // namespace dsrc

#endif
