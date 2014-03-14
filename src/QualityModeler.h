#ifndef QUALITYMODELER_H
#define QUALITYMODELER_H

#include "../include/dsrc/Globals.h"

#include "Common.h"
#include "RecordsProcessor.h"

namespace dsrc
{

namespace comp
{

class IQualityModeler
{
public:
	static const uchar HashSymbolNormal = IRecordsProcessor::HashSymbolNormal;
	static const uchar HashSymbolQuantized = IRecordsProcessor::HashSymbolQuantized;

	virtual ~IQualityModeler() {}

	virtual void ProcessStats(const QualityStats& stats_) = 0;
	virtual void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_) = 0;
	virtual void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_) = 0;

protected:
	static const uint32 MaxSymbolCount = 255;
};

} // namespace comp

} // namespace dsrc

#endif // QUALITYMODELER_H
