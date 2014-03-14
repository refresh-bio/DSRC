#ifndef DNAMODELHUFFMAN_H
#define DNAMODELHUFFMAN_H

#include "../include/dsrc/Globals.h"

#include "DnaModeler.h"
#include "Stats.h"
#include "BitMemory.h"
#include "Fastq.h"

#include "huffman.h"

namespace dsrc
{

namespace comp
{

class DnaModelerHuffman : public IDnaModeler
{
public:
	DnaModelerHuffman()
		:	symbolCount(0)
	{
		std::fill(symbols, symbols + MaxSymbolCount, 0);
	}

	void ProcessStats(const DnaStats &stats_);
	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_);
	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_);


private:
	static const uint32 MaxSymbolCount = DnaStats::MaxSymbolCount;

	uint32 symbolCount;
	uchar symbols[MaxSymbolCount];

	HuffmanEncoder coder;
};

} // namespace comp

} // namespace dsrc

#endif // DNAMODELHUFFMAN_H
