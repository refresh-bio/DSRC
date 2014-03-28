/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_QUALITYMODELERPROXY
#define H_QUALITYMODELERPROXY

#include "../include/dsrc/Globals.h"

#include "Fastq.h"
#include "QualityModeler.h"
#include "QualityPositionModeler.h"
#include "QualityRLEModeler.h"
#include "QualityOrderModeler.h"

namespace dsrc
{

namespace comp
{

class IQualityModelerProxy : public IQualityModeler
{
public:
	IQualityModelerProxy()
		:	modeler(NULL)
		,	currentSchemeId(SchemeNone)
	{}

	void ProcessStats(const QualityStats& stats_)
	{
		ASSERT(stats_.symbolCount < MaxSymbolCount);

		currentSchemeId = SelectSchemeId(stats_);
		if (currentSchemeId == SchemeNone)
			return;

		modeler = SelectModeler(currentSchemeId);
		ASSERT(modeler != NULL);
		modeler->ProcessStats(stats_);
	}

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_)
	{
		writer_.PutByte(currentSchemeId);

		if (currentSchemeId == SchemeNone)
			return;

		modeler = SelectModeler(currentSchemeId);
		ASSERT(modeler != NULL);
		modeler->Encode(writer_, records_, recordsCount_);
	}

	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_)
	{
		currentSchemeId = reader_.GetByte();

		if (currentSchemeId == SchemeNone)
			return;

		modeler = SelectModeler(currentSchemeId);
		ASSERT(modeler != NULL);
		modeler->Decode(reader_, records_, recordsCount_);
	}

protected:
	typedef byte SchemeId;

	static const uint32 MaxSymbolCount = 2566;
	static const SchemeId SchemeNone = 255;

	IQualityModeler* modeler;
	SchemeId currentSchemeId;

	virtual SchemeId SelectSchemeId(const QualityStats& stats_) = 0;
	virtual IQualityModeler* SelectModeler(SchemeId schemeId_) = 0;
};


class QualityNormalModelerProxy : public IQualityModelerProxy
{
public:
	QualityNormalModelerProxy(bool useQuantizedValues = false)
	{
		modelers[QualityPlain] = new QualityPositionModelerPlain(useQuantizedValues);
		modelers[QualityTruncated] = new QualityPositionModelerTruncated(useQuantizedValues);
		modelers[QualityRle] = new QualityRLEModeler(useQuantizedValues);
	}

	~QualityNormalModelerProxy()
	{
		delete modelers[QualityPlain];
		delete modelers[QualityTruncated];
		delete modelers[QualityRle];
	}

private:
	enum QualityLosslessSchemeEnum
	{
		QualityPlain = 0,
		QualityTruncated,
		QualityRle
	};

	IQualityModeler* modelers[3];

	SchemeId SelectSchemeId(const QualityStats& stats_)
	{
		if ((float)stats_.thLength / stats_.rleLength > 1.25f)
			return QualityRle;

		if ((float)stats_.rawLength / stats_.thLength > 1.10f)	// 1.10 is optimal, 1.25 results already in worse ratio
			return QualityTruncated;

		return QualityPlain;
	}

	IQualityModeler* SelectModeler(SchemeId schemeId_)
	{
		return modelers[schemeId_];
	}
};

class QualityOrderModelerProxyLossy : public IQualityModeler
{
public:
	QualityOrderModelerProxyLossy(uint32 order_)
		:	modeler(NULL)
	{
		ASSERT(order_ > 0 && order_ < 10);
		modeler = CreateModeler(order_);
	}

	~QualityOrderModelerProxyLossy()
	{
		delete modeler;
	}

	void ProcessStats(const QualityStats& stats_)
	{
		ASSERT(stats_.symbolCount < MaxSymbolCount);
		modeler->ProcessStats(stats_);
	}

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_)
	{
		modeler->Encode(writer_, records_, recordsCount_);
	}

	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_)
	{
		modeler->Decode(reader_, records_, recordsCount_);
	}

private:
	IQualityModeler* modeler;

	IQualityModeler* CreateModeler(uint32 order_)
	{
		switch (order_)
		{
			case 1:	return new TQualityLossyOrderPositionalModeler<8, 1>();
			case 2:	return new TQualityLossyOrderPositionalModeler<8, 2>();
			case 3:	return new TQualityLossyOrderPositionalModeler<8, 3>();
			case 4:	return new TQualityLossyOrderPositionalModeler<8, 4>();
			case 5:	return new TQualityLossyOrderPositionalModeler<8, 5>();
			case 6:	return new TQualityLossyOrderPositionalModeler<8, 6>();
			case 7:	return new TQualityLossyOrderPositionalModeler<8, 7>();
			case 8:	return new TQualityLossyOrderPositionalModeler<8, 8>();
			case 9:	return new TQualityLossyOrderPositionalModeler<8, 9>();
		}

		return NULL;
	}
};


class QualityOrderModelerProxyLossless : public IQualityModelerProxy
{
public:
	QualityOrderModelerProxyLossless(uint32 order_)
		:	order(order_)
	{
		ASSERT(order > 0 && order <= 2);

		std::fill(modelers, modelers + ModelersCount, (IQualityModeler*)0);
		//modeler = SelectModeler(QualitySym64);
	}

	~QualityOrderModelerProxyLossless()
	{
		for (uint32 i = 0; i < ModelersCount; ++i)
		{
			if (modelers[i] != NULL)
				delete modelers[i];
		}
	}

private:
	static const uint32 ModelersCount = 4 * 2;

	enum SymbolSchemeEnum
	{
		QualitySym16S = 0,
		QualitySym32S,
		QualitySym64S,
		QualitySym128S,

		QualitySym16F,
		QualitySym32F,
		QualitySym64F,
		QualitySym128F
	};

	const uint32 order;

	IQualityModeler* modelers[ModelersCount];

	IQualityModeler* CreateModeler(SchemeId scheme_)
	{
		if (order == 1)
		{
			switch (scheme_)
			{
				case QualitySym16S:		return new TQualityLosslessOrderTranslationalModeler<16, 3, 8>();
				case QualitySym32S:		return new TQualityLosslessOrderTranslationalModeler<32, 2, 8>();
				case QualitySym64S:		return new TQualityLosslessOrderTranslationalModeler<64, 1, 8>();
				case QualitySym128S:	return new TQualityLosslessOrderTranslationalModeler<128, 1, 8>();

				case QualitySym16F:		return new TQualityLosslessOrderTranslationalModeler<16, 3, 16>();
				case QualitySym32F:		return new TQualityLosslessOrderTranslationalModeler<32, 2, 32>();
				case QualitySym64F:		return new TQualityLosslessOrderTranslationalModeler<64, 1, 64>();
				case QualitySym128F:	return new TQualityLosslessOrderTranslationalModeler<128, 1, 128>();
			}
		}
		else
		{
			switch (scheme_)
			{
				case QualitySym16S:		return new TQualityLosslessOrderTranslationalModeler<16, 4, 8>();		// can try one up
				case QualitySym32S:		return new TQualityLosslessOrderTranslationalModeler<32, 3, 8>();
				case QualitySym64S:		return new TQualityLosslessOrderTranslationalModeler<64, 2, 8>();
				case QualitySym128S:	return new TQualityLosslessOrderTranslationalModeler<128, 1, 8>();		// can try one up

				case QualitySym16F:		return new TQualityLosslessOrderTranslationalModeler<16, 4, 16>();		// can try one up
				case QualitySym32F:		return new TQualityLosslessOrderTranslationalModeler<32, 3, 32>();
				case QualitySym64F:		return new TQualityLosslessOrderTranslationalModeler<64, 2, 64>();
				case QualitySym128F:	return new TQualityLosslessOrderTranslationalModeler<128, 1, 128>();	// can try one up
			}
		}

		return NULL;
	}

	SchemeId SelectSchemeId(const QualityStats &stats_)
	{
		SchemeId scheme = SchemeNone;
		for (uint32 i = 0; i < ModelersCount; ++i)
		{
			if ((16U << i) >= stats_.symbolCount)
			{
				scheme = i;
				break;
			}
		}

		// use special scheme modification for const-length records with RLE factor > 1.175 and in order 2
		if (scheme != SchemeNone && order == 2)
		{
			double rleRatio = (double)stats_.rawLength / stats_.rleLength;
			if (stats_.maxLength == stats_.minLength && rleRatio > 1.175)
				scheme += 4;	// offset
		}

		return scheme;
	}

	IQualityModeler* SelectModeler(SchemeId scheme_)
	{
		if (scheme_ == SchemeNone)
			return NULL;

		if (modelers[scheme_] == NULL)
			modelers[scheme_] = CreateModeler(scheme_);
		return modelers[scheme_];
	}
};

} // namespace comp

} // namespace dsrc

#endif // H_QUALITYMODELERPROXY
