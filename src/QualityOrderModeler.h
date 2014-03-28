/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
#ifndef H_QUALITYORDERMODELER
#define H_QUALITYORDERMODELER

#include "../include/dsrc/Globals.h"

#include "QualityEncoder.h"
#include "Fastq.h"
#include "QualityModeler.h"

#include "huffman.h"
#include "utils.h"

namespace dsrc
{

namespace comp
{

template <class _TQualityEncoder>
class TQualityOrderModeler : public IQualityModeler
{
public:
	void ProcessStats(const QualityStats& stats_)
	{
		encoder.ProcessStats(stats_);
	}

	void Encode(core::BitMemoryWriter& writer_, const fq::FastqRecord* records_, uint32 recordsCount_)
	{
		encoder.Store(writer_);

		model.Clear();

		RangeEncoder coder(writer_);
		coder.Start();

		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			const fq::FastqRecord& r = records_[i];
			encoder.Encode(r, model, coder);
		}
		coder.End();
	}

	void Decode(core::BitMemoryReader& reader_, fq::FastqRecord* records_, uint32 recordsCount_)
	{
		encoder.Read(reader_);

		model.Clear();

		RangeDecoder coder(reader_);
		coder.Start();

		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			fq::FastqRecord& r = records_[i];
			encoder.Decode(r, model, coder);
		}
		coder.End();
	}

private:
	typedef _TQualityEncoder Encoder;
	typedef typename Encoder::Model Model;

	Model model;
	Encoder encoder;
};



template <uint32 _TSymbolCount, uint32 _TOrder>
class TQualityLossyOrderNormalModeler : public TQualityOrderModeler<
													TNormalQualityEncoder<
														TQualityModel<_TSymbolCount, _TOrder>,
														SpecialSymbolLossyHandlerStrategy>
												>
{
private:
	typedef TQualityOrderModeler<TNormalQualityEncoder<
									TQualityModel<_TSymbolCount, _TOrder>,
									SpecialSymbolLossyHandlerStrategy>
								> Super;

public:
	using Super::Encode;
	using Super::Decode;
	using Super::ProcessStats;

};

template <uint32 _TSymbolCount, uint32 _TOrder>
class TQualityLossyOrderPositionalModeler : public TQualityOrderModeler<
														TPositionalQualityEncoder<
															TQualityModelExt<_TSymbolCount, _TOrder>,
															SpecialSymbolLossyHandlerStrategy>
													>
{
private:
	typedef TQualityOrderModeler<TPositionalQualityEncoder<
									TQualityModelExt<_TSymbolCount, _TOrder>,
									SpecialSymbolLossyHandlerStrategy>
								> Super;

public:
	using Super::Encode;
	using Super::Decode;
	using Super::ProcessStats;
};

/*
template <uint32 _TSymbolCount, uint32 _TOrder>
class TQualityLosslessOrderNormalModeler : public TQualityOrderModeler<
														TNormalQualityEncoder<
															TQualityModel<_TSymbolCount, _TOrder>,
															SpecialSymbolNormalHandlerStrategy>
													>
{
private:
	typedef TQualityOrderModeler<TNormalQualityEncoder<
									TQualityModel<_TSymbolCount, _TOrder>,
									SpecialSymbolNormalHandlerStrategy>
								> Super;

public:
	using Super::Encode;
	using Super::Decode;
	using Super::ProcessStats;
};
*/

template <uint32 _TSymbolCount, uint32 _TOrder, uint32 _TSymbolRescale = 8>
class TQualityLosslessOrderTranslationalModeler : public TQualityOrderModeler<
															TTranslationalQualityEncoder<
																TQualityModelExt<_TSymbolCount, _TOrder>,
																SpecialSymbolNormalHandlerStrategy,
																_TSymbolRescale >
														>
{
private:
	typedef TQualityOrderModeler<TTranslationalQualityEncoder<
									TQualityModelExt<_TSymbolCount, _TOrder>,
									SpecialSymbolNormalHandlerStrategy,
									_TSymbolRescale >
								> Super;

public:
	using Super::Encode;
	using Super::Decode;
	using Super::ProcessStats;
};

} // namespace comp

} // namespace dsrc

#endif // H_QUALITYORDERMODELER
