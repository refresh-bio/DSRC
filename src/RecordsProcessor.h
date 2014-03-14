/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_RECORDSPROCESSOR
#define H_RECORDSPROCESSOR

#include "../include/dsrc/Globals.h"

#include <algorithm>

#include "Fastq.h"
#include "Stats.h"
#include "Crc32.h"

namespace dsrc
{

namespace comp
{


class FastqChecksumHasher
{
public:
	FastqChecksumHasher()
		:	flags(fq::FastqChecksum::CALC_NONE)
	{}

	void SetFlags(uint32& flags_)
	{
		flags = flags_;
	}

	void Update(const fq::FastqRecord& r_)
	{
		if (flags & fq::FastqChecksum::CALC_TAG)
			hashers[0].UpdateCrc(r_.title, r_.titleLen);
		if (flags & fq::FastqChecksum::CALC_SEQUENCE)
			hashers[1].UpdateCrc(r_.sequence, r_.sequenceLen);
		if (flags & fq::FastqChecksum::CALC_QUALITY)
			hashers[2].UpdateCrc(r_.quality, r_.qualityLen);
	}

	void Reset()
	{
		for (uint32 i = 0; i < 3; ++i)
			hashers[i].Reset();
	}

	fq::FastqChecksum GetChecksum() const
	{
		fq::FastqChecksum checksum;
		checksum.tag = hashers[0].GetHash();
		checksum.sequence = hashers[1].GetHash();
		checksum.quality = hashers[2].GetHash();
		return checksum;
	}

private:
	uint32 flags;
	core::Crc32Hasher hashers[3];
};


class IRecordsProcessor
{
public:
	typedef fq::FastqChecksum::CheksumFlags ChecksumFlags;

	static const uchar HashSymbolNormal = 2;		// not '0'!!!
	static const uchar HashSymbolQuantized = 1;

	IRecordsProcessor(uint32 qualityOffset_ = 33, bool colorSpace_ = false)
		:	qualityOffset(qualityOffset_)
		,	colorSpace(colorSpace_)
	{
		ASSERT(qualityOffset_ >= 33 && qualityOffset_ <= 64);
		// offset = 33 - standard Sanger + new Illumina
		// offset = 59 - SOLEXA
		// offset = 64 - old Illumina
	}

	virtual ~IRecordsProcessor()
	{}

	void ProcessFromColorSpace(fq::FastqRecord& rec_)
	{
		ProcessRecordFromColorSpace(rec_);

		// update stats
		//
		if (csStats.seqBegin == ColorSpaceStats::EmptySymbol)
		{
			csStats.seqBegin = rec_.sequence[0];
			csStats.quaBegin = rec_.quality[0];
		}
		csStats.constBeginSym &= csStats.seqBegin == rec_.sequence[0];
		ASSERT(!csStats.constBeginSym || csStats.quaBegin == rec_.quality[0]);
	}

	void ProcessToColorSpace(fq::FastqRecord& rec_, uchar seq0, uchar qua0)
	{
		ProcessRecordToColorSpace(rec_, csStats.constBeginSym, seq0, qua0);
	}

	virtual void InitializeStats();
	virtual void FinalizeStats();

	// returns CRC32 chekcsum if enabled
	fq::FastqChecksum ProcessForward(fq::FastqRecord* records_, uint64 recordsCount_, uint32 flags_ = fq::FastqChecksum::CALC_NONE);
	fq::FastqChecksum ProcessBackward(fq::FastqRecord* records_, uint64 recordsCount_, uint32 flags_ = fq::FastqChecksum::CALC_NONE);

	const DnaStats& GetDnaStats() const
	{
		return dnaStats;
	}

	const QualityStats& GetQualityStats() const
	{
		return qualityStats;
	}

	const ColorSpaceStats& GetColorSpaceStats() const
	{
		return csStats;
	}

	void SetColorSpaceStats(const ColorSpaceStats& stats_)
	{
		csStats = stats_;
	}


protected:
	static const uchar InvalidValue = 255;

	const uint32 qualityOffset;
	const bool colorSpace;

	DnaStats dnaStats;
	QualityStats qualityStats;
	ColorSpaceStats csStats;

	FastqChecksumHasher fastqHasher;

private:
	void ProcessRecordFromColorSpace(fq::FastqRecord& rec_);
	void ProcessRecordToColorSpace(fq::FastqRecord& rec_, bool useConstDelta_, uchar seqStart_, uchar quaStart_);

	virtual void ProcessForward(fq::FastqRecord &rec_) = 0;
	virtual void ProcessBackward(fq::FastqRecord &rec_) = 0;
};

class LosslessRecordsProcessor : public IRecordsProcessor
{
public:
	LosslessRecordsProcessor(uint32 qualityOffset_, bool colorSpace_);

protected:
	uchar dnaToIndexTable[128];
	uchar dnaFromIndexTable[DnaStats::MaxSymbolCount];

private:
	void ProcessForward(fq::FastqRecord &rec_);
	void ProcessBackward(fq::FastqRecord &rec_);
};


class LossyRecordsProcessor : public LosslessRecordsProcessor
{
public:
	LossyRecordsProcessor(uint32 qualityOffset_, bool colorSpace_);

	void FinalizeStats()
	{
		IRecordsProcessor::FinalizeStats();

		ASSERT(dnaStats.symbolCount <= 4);
		ASSERT(qualityStats.symbolCount <= 8);
	}

private:
	uchar qualityToIndexTable[64];
	uchar qualityFromIndexTable[8];

	void ProcessForward(fq::FastqRecord &rec_);
	void ProcessBackward(fq::FastqRecord &rec_);
};

} // namespace comp

} // namespace dsrc


#endif // H_RECORDSPROCESSOR
