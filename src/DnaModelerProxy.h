/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DNAMODELERPROXY
#define H_DNAMODELERPROXY

#include "../include/dsrc/Globals.h"

#include "Stats.h"
#include "DnaModelerBasicB2.h"
#include "DnaModelerHuffman.h"
#include "DnaModelerRCO.h"

#include <algorithm>


namespace dsrc
{

namespace comp
{

class IDnaModelerProxy : public IDnaModeler
{
public:
	IDnaModelerProxy()
		:	modeler(NULL)
		,	currentSchemeId(SchemeNone)
	{}

	void ProcessStats(const DnaStats& stats_)
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

	static const uint32 MaxSymbolCount = 20;
	static const SchemeId SchemeNone = 255;

	IDnaModeler* modeler;
	SchemeId currentSchemeId;

	virtual SchemeId SelectSchemeId(const DnaStats& stats_) = 0;
	virtual IDnaModeler* SelectModeler(SchemeId schemeId_) = 0;
};


class DnaNormalModelerProxy : public IDnaModelerProxy
{
private:
	static const uint32 MaxSymbolCount = DnaStats::MaxSymbolCount;

	enum Order0Schemes
	{
		SchemeB2 = 0,
		SchemeHuffman
	};

	DnaModelerBasicB2 basicModeler;
	DnaModelerHuffman huffModeler;

	SchemeId SelectSchemeId(const DnaStats &stats_)
	{
		if (stats_.symbolCount == 0)
			return SchemeNone;

		if (stats_.symbolCount <= 4)
			return SchemeB2;

		return SchemeHuffman;
	}

	IDnaModeler* SelectModeler(SchemeId schemeId_)
	{
		if (schemeId_ == SchemeB2)
			return &basicModeler;

		if (schemeId_ == SchemeHuffman)
			return &huffModeler;

		return NULL;
	}
};

class DnaOrderModelerProxy : public IDnaModelerProxy
{
public:
	DnaOrderModelerProxy(uint32 order_)
		:	order(order_)
		,	modeler4s(NULL)
		,	modeler8s(NULL)
	{
		ASSERT(order_ > 0 && order_ < 10);

		 modeler4s = CreateModeler(Scheme4Sym);
	}

	~DnaOrderModelerProxy()
	{
		if (modeler4s != NULL)
			delete modeler4s;

		if (modeler8s != NULL)
			delete modeler8s;
	}

private:
	static const uint32 MaxSymbolCount = 8;
	const uint32 order;

	enum OrderNSchemes
	{
		Scheme4Sym = 0,
		Scheme8Sym
	};

	IDnaModeler* modeler4s;
	IDnaModeler* modeler8s;

	SchemeId SelectSchemeId(const DnaStats &stats_)
	{
		if (stats_.symbolCount == 0)
			return SchemeNone;

		if (stats_.symbolCount <= 4)
			return Scheme4Sym;

		ASSERT(stats_.symbolCount <= 8);
		return Scheme8Sym;
	}

	IDnaModeler* SelectModeler(SchemeId scheme_)
	{
		if (scheme_ == Scheme4Sym)
		{
			if (modeler4s == NULL)
				modeler4s = CreateModeler(Scheme4Sym);
			return modeler4s;
		}

		if (scheme_ == Scheme8Sym)
		{
			if (modeler8s == NULL)
				modeler8s = CreateModeler(Scheme8Sym);
			return modeler8s;
		}

		ASSERT(0);
		return NULL;
	}

	IDnaModeler* CreateModeler(SchemeId scheme_)
	{
		if (scheme_ == Scheme4Sym)
		{
			switch (order)
			{
				case 1: return new TDnaRCOrderModeler<1, 4>();
				case 2: return new TDnaRCOrderModeler<2, 4>();
				case 3: return new TDnaRCOrderModeler<3, 4>();
				case 4: return new TDnaRCOrderModeler<4, 4>();
				case 5: return new TDnaRCOrderModeler<5, 4>();
				case 6: return new TDnaRCOrderModeler<6, 4>();
				case 7: return new TDnaRCOrderModeler<7, 4>();
				case 8: return new TDnaRCOrderModeler<8, 4>();
				case 9: return new TDnaRCOrderModeler<9, 4>();
			}
		}

		if (scheme_ == Scheme8Sym)
		{
			switch (order)
			{
				case 1: return new TDnaRCOrderModeler<1, 8>();
				case 2: return new TDnaRCOrderModeler<2, 8>();
				case 3: return new TDnaRCOrderModeler<3, 8>();
				case 4: return new TDnaRCOrderModeler<4, 8>();
				case 5: return new TDnaRCOrderModeler<5, 8>();
				case 6: return new TDnaRCOrderModeler<6, 8>();
				case 7: return new TDnaRCOrderModeler<7, 8>();	// lock on the 7th order due to too high memory usage
				case 8: return new TDnaRCOrderModeler<7, 8>();	// on higher orders : 7th - 2^25, 8th -2^28 ...
				case 9: return new TDnaRCOrderModeler<7, 8>();
			}
		}

		ASSERT(0);
		return NULL;
	}
};

} // namespace comp

} // namespace dsrc

#endif // H_DNAMODELERPROXY
