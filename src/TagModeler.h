/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
  
  Version: 2.00
*/

#ifndef H_TAGMODELER
#define H_TAGMODELER

#include "../include/Globals.h"

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


class TagModeler
{
public:
	TagModeler()
		:	recordCounter(0)
	{}

	void InitializeFieldsStats(const fq::FastqRecord& rec_);
	void UpdateFieldsStats(const fq::FastqRecord& rec_);
	void FinalizeFieldsStats();

	void StartEncoding(core::BitMemoryWriter& writer_);
	void EncodeNextFields(core::BitMemoryWriter& writer_, const fq::FastqRecord& rec_);
	void FinishEncoding(core::BitMemoryWriter& writer_);

	void StartDecoding(core::BitMemoryReader& reder_);
	void DecodeNextFields(core::BitMemoryReader& reader_, fq::FastqRecord& rec_);
	void FinishDecoding(core::BitMemoryReader& reader_);

private:
	static const uint32 MAX_FIELD_STAT_LEN = 128;
	static const uint32 MAX_NUM_VAL_HUF	= 512;

	std::vector<Field> fields;
	std::vector<int32> prevFieldValues;
	uint32 recordCounter;

	void UpdateNumericField(Field& field_, int32 curValue_, int32 prevValue_);

	void StoreFields(core::BitMemoryWriter &bit_stream);
	void StoreNumericField(core::BitMemoryWriter& bit_memory, Field& field_, int32 curValue_, int32 prevValue_);

	void ReadFields(core::BitMemoryReader &bit_stream);
	uint32 ReadNumericField(core::BitMemoryReader& bit_stream, Field& field_, int32 prevValue_);
};

} // namespace comp

} // namespace dsrc

#endif // TAGMODEL_H
