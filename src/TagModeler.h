/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_TAGMODELER
#define H_TAGMODELER

#include "../include/dsrc/Globals.h"

#include <vector>
#include <map>

#include "Common.h"

#include "huffman.h"
#include "utils.h"

namespace dsrc
{

namespace comp
{

// a lot of algorithms code in this module was taken from DSRC1.0 and needs a-bit-o'-refactoring

struct Field
{
	static const uint32 HUF_GLOBAL_SIZE = 512;
	static const uint32 HUF_LOCAL_SIZE = 256;
	static const uint32 MAX_FIELD_STAT_LEN = 128;
	static const uint32 MAX_NUM_VAL_HUF	= 512;

	uint32 len;
	uint32 min_len;
	uint32 max_len;
	uchar sep;
	bool is_constant;
	bool is_len_constant;
	bool is_numeric;
	int32 min_value;
	int32 max_value;
	int32 min_delta;
	int32 max_delta;
	uint32 no_of_bits_per_num;
	uint32 no_of_bits_per_value;
	uint32 no_of_bits_per_len;
	bool is_delta_coding;

	struct RleInfo
	{
		int32		cur_sym;
		uint32		cur_len;
		uint32		run_len;
		std::vector<uchar> lens;
		RleInfo()
			:	cur_sym(0)
			,	cur_len(0)
			,	run_len(0)
		{}
	};
	RleInfo rle_val;		// todo: bool flags to control collecting stats
	RleInfo rle_delta;
	bool try_rle_val;
	bool try_rle_delta;
	bool is_delta_const;
	bool var_stat_encode;

	enum NumericSchemeEnum { None = 0, ValueVar, ValueRle, DeltaVar, DeltaRle, DeltaConst };
	uchar numeric_scheme;

	uchar *data;
	uint32 *global_stats;
	uint32 **stats;
	uint32 *raw_stats;
	bool *ham_mask;

	HuffmanEncoder *Huffman_global;
	std::vector<HuffmanEncoder*> Huffman_local;

	std::map<int32, int32> num_values;
	std::map<int32, int32> delta_values;

	std::vector<std::map<char,uint32> > chars;	// @ pos [char, count]

	Field();
	Field(const Field& f_);

	~Field();

	bool IsNum(uint32& num_) const { return core::is_num(data, len, num_); }
	//uint32 ToNum() const { return utils::to_num(data, len); }
	uint32 ToString() const { return core::to_string(data, len); }
};


struct TagStats
{
	static const uint32 MaxSymbolCount = 128;

	std::vector<Field> fields;
	uint32 minTitleLen;
	uint32 maxTitleLen;
	uint32 symbolFreqs[MaxSymbolCount];
	bool mixedFormatting;

	TagStats()
	{
		Reset();
	}

	void Reset()
	{
		fields.clear();				// TODO: reset fields, do not delete!
		minTitleLen = 0xFFFFFFFF;
		maxTitleLen = 0;
		std::fill(symbolFreqs, symbolFreqs + MaxSymbolCount, 0);
		mixedFormatting = false;
	}
};


class TagAnalyzer
{
public:
	void InitializeFieldsStats(const fq::FastqRecord& rec_);
	void UpdateFieldsStats(const fq::FastqRecord& rec_);
	void FinalizeFieldsStats();

	TagStats& GetStats()
	{
		return stats;
	}

private:
	TagStats stats;
	std::vector<int32> prevFieldValues;
	uint32 recordCounter;

	void UpdateNumericField(Field& field_, int32 curValue_, int32 prevValue_);
};


class ITagEncoder
{
public:
	ITagEncoder()
	{}

	virtual ~ITagEncoder() {}

	virtual void StartEncoding(core::BitMemoryWriter& writer_, TagStats* stats_) = 0;
	virtual void EncodeNextFields(core::BitMemoryWriter& writer_, const fq::FastqRecord& rec_) = 0;
	virtual void FinishEncoding(core::BitMemoryWriter& writer_) = 0;
};


class ITagDecoder
{
public:
	virtual ~ITagDecoder() {}

	virtual void StartDecoding(core::BitMemoryReader& reader_) = 0;
	virtual void DecodeNextFields(core::BitMemoryReader& reader_, fq::FastqRecord& rec_) = 0;
	virtual void FinishDecoding(core::BitMemoryReader& reader_) = 0;
};


class ITagTokenizer
{
public:
	ITagTokenizer()
		:	recordCounter(0)
	{}

protected:
	std::vector<int32> prevFieldValues;
	uint32 recordCounter;
};


class TagTokenizerEncoder : public ITagTokenizer, public ITagEncoder
{
public:
	TagTokenizerEncoder()
		:	stats(NULL)
	{}

	void StartEncoding(core::BitMemoryWriter& writer_, TagStats* stats_);
	void EncodeNextFields(core::BitMemoryWriter& writer_, const fq::FastqRecord& rec_);
	void FinishEncoding(core::BitMemoryWriter& writer_);

private:
	TagStats* stats;

	void StoreFields(core::BitMemoryWriter &bit_stream);
	void StoreNumericField(core::BitMemoryWriter& bit_memory, Field& field_, int32 curValue_, int32 prevValue_);
};


class TagTokenizerDecoder : public ITagTokenizer, public ITagDecoder
{
public:
	void StartDecoding(core::BitMemoryReader& reader_);
	void DecodeNextFields(core::BitMemoryReader& reader_, fq::FastqRecord& rec_);
	void FinishDecoding(core::BitMemoryReader& reader_);

private:
	TagStats stats;

	void ReadFields(core::BitMemoryReader &bit_stream);
	uint32 ReadNumericField(core::BitMemoryReader& bit_stream, Field& field_, int32 prevValue_);
};


class ITagRawCoder
{
public:
	ITagRawCoder()
		:	symbolCount(0)
		,	titleLenBits(0)
		,	encoder(NULL)
	{}

	virtual ~ITagRawCoder()
	{
		if (encoder)
			delete encoder;
	}

protected:
	static const uint32 MaxSymbolCount = 128;
	static const char EmptySymbol = -1;

	uint32 symbolCount;
	char symbols[MaxSymbolCount];
	uint32 titleLenBits;
	HuffmanEncoder* encoder;		// try non-dynamic + Reset()
};


class TagRawEncoder : public ITagEncoder, ITagRawCoder
{
public:
	TagRawEncoder()
		:	stats(NULL)
	{}

	void StartEncoding(core::BitMemoryWriter& writer_, TagStats* stats_);
	void EncodeNextFields(core::BitMemoryWriter& writer_, const fq::FastqRecord& rec_);
	void FinishEncoding(core::BitMemoryWriter& writer_);

private:
	TagStats* stats;
};


class TagRawDecoder : public ITagDecoder, ITagRawCoder
{
public:
	TagRawDecoder()
		:	minTitleLen(0)
		,	maxTitleLen(0)
	{}

	void StartDecoding(core::BitMemoryReader& reader_);
	void DecodeNextFields(core::BitMemoryReader& reader_, fq::FastqRecord& rec_);
	void FinishDecoding(core::BitMemoryReader& reader_);

private:
	uint32 minTitleLen;
	uint32 maxTitleLen;
};


class TagModeler
{
public:
	enum TagEncodingScheme
	{
		TagTokenizeHuffman = 0,
		TagRawHuffman
	};

	TagModeler()
		:	analyzer(NULL)
	{
		encoders[0] = encoders[1] = NULL;
		decoders[0] = decoders[1] = NULL;
	}

	~TagModeler()
	{
		if (encoders[0] != NULL) delete encoders[0];
		if (encoders[1] != NULL) delete encoders[1];
		if (decoders[0] != NULL) delete decoders[0];
		if (decoders[1] != NULL) delete decoders[1];
	}

	ITagEncoder* SelectEncoder(TagEncodingScheme scheme_)
	{
		if (encoders[scheme_] == NULL)
		{
			if (scheme_ == TagTokenizeHuffman)
				encoders[scheme_] = new TagTokenizerEncoder();
			else
				encoders[scheme_] = new TagRawEncoder();
		}
		return encoders[scheme_];
	}

	ITagDecoder* SelectDecoder(TagEncodingScheme scheme_)
	{
		if (decoders[scheme_] == NULL)
		{
			if (scheme_ == TagTokenizeHuffman)
				decoders[scheme_] = new TagTokenizerDecoder();
			else
				decoders[scheme_] = new TagRawDecoder();
		}
		return decoders[scheme_];
	}

	TagAnalyzer* GetAnalyzer()
	{
		if (analyzer == NULL)
			analyzer = new TagAnalyzer();
		return analyzer;
	}

private:
	TagAnalyzer* analyzer;
	ITagEncoder* encoders[2];
	ITagDecoder* decoders[2];
};

} // namespace comp

} // namespace dsrc

#endif // TAGMODEL_H
