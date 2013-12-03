/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "DnaModelerHuffman.h"

namespace dsrc
{

namespace comp
{

using namespace core;
using namespace fq;

void DnaModelerHuffman::ProcessStats(const DnaStats &stats_)
{
	// copy symbol info
	//
	ASSERT(stats_.symbolCount > 0);
	symbolCount = stats_.symbolCount;
	std::copy(stats_.symbols, stats_.symbols + MaxSymbolCount, symbols);


	// prepare Huffman
	//
	coder.Restart(symbolCount);

	for (uint32 i = 0; i < symbolCount; ++i)
	{
		coder.Insert(stats_.symbolFreqs[symbols[i]]);
	}
	coder.Complete();
}

void DnaModelerHuffman::Encode(BitMemoryWriter& writer_, const FastqRecord* records_, uint32 recordsCount_)
{
	ASSERT(symbolCount > 0);

	// save codes
	//
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		writer_.PutBit(symbols[i] != DnaStats::EmptySymbol);
	}
	writer_.FlushPartialWordBuffer();

	// save tree
	//
	coder.StoreTree(writer_);

	// start encoding
	//
	const HuffmanEncoder::Code* symbolCodes = coder.GetCodes();

	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		const FastqRecord& r = records_[i];

		for (uint32 j = 0; j < r.sequenceLen; ++j)
		{
			const HuffmanEncoder::Code& code = symbolCodes[symbols[r.sequence[j]]];
			writer_.PutBits(code.code, code.len);
		}
	}

	writer_.FlushPartialWordBuffer();
}

void DnaModelerHuffman::Decode(BitMemoryReader& reader_, FastqRecord* records_, uint32 recordsCount_)
{
	// load sybols
	//
	symbolCount = 0;
	std::fill(symbols, symbols + MaxSymbolCount, +DnaStats::EmptySymbol);
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		if (reader_.GetBit())
		{
			symbols[symbolCount++] = i;
		}
	}

	ASSERT(symbolCount > 0);

	// load huffman
	//
	coder.LoadTree(reader_);

	for (uint32 i = 0; i < recordsCount_; ++i)
	{
		FastqRecord& r = records_[i];

		for (uint32 j = 0; j < r.sequenceLen; ++j)
		{
			uint32 bit = reader_.GetBits(coder.GetMinLen());
			int32 sidx = coder.DecodeFast(bit);
			while (sidx < 0)
			{
				bit = reader_.GetBit();
				sidx = coder.Decode(bit);
			};

			r.sequence[j] = symbols[sidx];
		}
	}

	reader_.FlushInputWordBuffer();
}

} // namespace comp

} // namespace dsrc
