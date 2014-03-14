/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_STATS
#define H_STATS

#include "../include/dsrc/Globals.h"

#include <algorithm>

namespace dsrc
{

namespace comp
{

struct ColorSpaceStats
{
	static const uint32 EmptySymbol = 255;

	bool constBeginSym;
	uchar seqBegin;
	uchar quaBegin;

	ColorSpaceStats()
	{
		Clear();
	}

	void Clear()
	{
		constBeginSym = true;
		seqBegin = EmptySymbol;
		quaBegin = EmptySymbol;
	}
};

struct DnaStats
{
	static const uint32 MaxSymbolCount = 20;
	static const uint32 EmptySymbol = 255;

	uint32 symbolCount;
	uint32 symbolFreqs[MaxSymbolCount];
	uchar symbols[MaxSymbolCount];

	DnaStats()
	{
		Clear();
	}

	void Clear()
	{
		symbolCount = 0;
		std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
		std::fill(symbolFreqs, symbolFreqs + MaxSymbolCount, 0);
	}
};

struct QualityStats
{
	static const uint32 MaxSymbolCount = 256;
	static const uint32 EmptySymbol = 255;

	uint32 symbolCount;
	uint32 symbolFreqs[MaxSymbolCount];
	uchar symbols[MaxSymbolCount];

	uint32 minLength;
	uint32 maxLength;

	uint32 rawLength;
	uint32 thLength;
	uint32 rleLength;

	uint32 symbolThreshold;

	QualityStats()
	{
		Clear();
	}

	void Clear()
	{
		symbolCount = 0;
		std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
		std::fill(symbolFreqs, symbolFreqs + MaxSymbolCount, 0);

		minLength = (uint32)-1;
		maxLength = 0;
		rawLength = thLength = rleLength = 0;

		symbolThreshold = 0;
	}
};

} // namespace comp

} // namespace dsrc

#endif // H_STATS
