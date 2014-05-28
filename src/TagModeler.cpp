/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "TagModeler.h"
#include "Fastq.h"

namespace dsrc
{

namespace comp
{

using namespace core;
using namespace fq;

// ********************************************************************************************
Field::Field()
	:	len(0)
	,	min_len(0)
	,	max_len(0)
	,	sep(0)
	,	is_constant(false)
	,	is_len_constant(false)
	,	is_numeric(false)
	,	min_value(1 << 30)
	,	max_value(-(1 << 30))
	,	min_delta(1 << 30)
	,	max_delta(-(1 << 30))
	,	no_of_bits_per_num(0)
	,	no_of_bits_per_value(0)
	,	no_of_bits_per_len(0)
	,	is_delta_coding(false)
	,	try_rle_val(false)
	,	try_rle_delta(false)
	,	is_delta_const(false)
	,	var_stat_encode(false)
	,	numeric_scheme(None)
	,	data(NULL)
	,	global_stats(NULL)
	,	stats(NULL)
	,	raw_stats(NULL)
	,	ham_mask(NULL)
	,	Huffman_global(NULL)
{}


// ********************************************************************************************
Field::Field(const Field &y)
	:	len(y.len)
	,	min_len(y.min_len)
	,	max_len(y.max_len)
	,	sep(y.sep)
	,	is_constant(y.is_constant)
	,	is_len_constant(y.is_len_constant)
	,	is_numeric(y.is_numeric)
	,	min_value(y.min_value)
	,	max_value(y.max_value)
	,	min_delta(y.min_delta)
	,	max_delta(y.max_delta)
	,	no_of_bits_per_num(y.no_of_bits_per_num)
	,	no_of_bits_per_value(y.no_of_bits_per_value)
	,	no_of_bits_per_len(y.no_of_bits_per_len)
	,	is_delta_coding(y.is_delta_coding)
	,	try_rle_val(y.try_rle_val)
	,	try_rle_delta(y.try_rle_delta)
	,	is_delta_const(y.is_delta_const)
	,	var_stat_encode(y.var_stat_encode)
	,	numeric_scheme(y.numeric_scheme)
{
	Field &x = const_cast<Field &>(y);

	if (x.data && len)
	{
		data   = x.data;
		x.data = NULL;
	}
	else
	{
		data = NULL;
	}

	if (x.ham_mask && len)
	{
		ham_mask   = x.ham_mask;
		x.ham_mask = NULL;
	}
	else
	{
		ham_mask = NULL;
	}

	if (x.stats)
	{
		stats       = x.stats;
		raw_stats   = x.raw_stats;
		x.stats     = NULL;
		x.raw_stats = NULL;
	}
	else
	{
		stats     = NULL;
		raw_stats = NULL;
	}

	if (x.global_stats)
	{
		global_stats = x.global_stats;
		x.global_stats = NULL;
	}
	else
	{
		global_stats = NULL;
	}

	if (x.Huffman_global)
	{
		Huffman_global = x.Huffman_global;
		x.Huffman_global = new HuffmanEncoder(Field::HUF_GLOBAL_SIZE);
	}
	else
	{
		Huffman_global = NULL;
	}
}


// ********************************************************************************************
Field::~Field()
{
	if (data)
		delete[] data;
	if (ham_mask)
		delete[] ham_mask;

	if (stats)
		delete[] stats;
	if (global_stats)
		delete[] global_stats;
	if (raw_stats)
		delete[] raw_stats;

	if (Huffman_global)
		delete Huffman_global;

	for (uint32 i = 0; i < Huffman_local.size(); ++i)
	{
		if (Huffman_local[i])
			delete Huffman_local[i];
	}
	Huffman_local.resize(0);
}

void TagAnalyzer::InitializeFieldsStats(const FastqRecord &rec_)
{
	const char *c_separators = " ._,=:/-#"; //9
	const std::vector<uchar> separators(c_separators, c_separators + 9 + 1);

	stats.Reset();
	prevFieldValues.clear();

	uint32 start_pos = 0;
	uint32 n_field = 0;

	for (uint32 i = 0; i <= rec_.titleLen; ++i)
	{
		stats.symbolFreqs[rec_.title[i]] += (i != rec_.titleLen);

		if (!std::count(separators.begin(), separators.end(), rec_.title[i]) && (i != rec_.titleLen))
			continue;

		stats.fields.push_back(Field());
		Field& field = stats.fields[n_field];

		field.data = new uchar[i-start_pos+1];
		std::copy(rec_.title + start_pos, rec_.title + i, field.data);

		field.data[i-start_pos] = '\0';
		field.len = i - start_pos;
		field.max_len = field.len;
		field.min_len = field.len;
		field.sep = rec_.title[i];
		field.is_constant = true;
		field.is_len_constant = true;

		uint32 num_val;
		field.is_numeric = field.IsNum(num_val);

		ASSERT(field.ham_mask == NULL);
		field.ham_mask = new bool[field.len];
		field.is_delta_const = false;
		field.try_rle_val = false;
		field.try_rle_delta	= false;
		field.var_stat_encode = false;
		field.numeric_scheme = Field::None;

		if (field.is_numeric)
		{
			field.min_value = num_val;
			field.max_value = field.min_value;
			field.num_values.clear();
			field.num_values[field.min_value]++;

			field.min_delta = (1 << 30);
			field.max_delta = -(1 << 30);
			field.delta_values.clear();
		}

		std::fill(field.ham_mask, field.ham_mask + field.len, true);

		start_pos = i + 1;
		n_field++;
	}

	// prepare counter for updating fields statistics
	//
	recordCounter = 0;
	prevFieldValues.resize(stats.fields.size(), 0);
}

void TagAnalyzer::UpdateFieldsStats(const FastqRecord &rec_)
{
	ASSERT(rec_.titleLen > 0);

	stats.minTitleLen = MIN(stats.minTitleLen, rec_.titleLen);
	stats.maxTitleLen = MAX(stats.maxTitleLen, rec_.titleLen);

	if (stats.mixedFormatting)
	{
		for (uint32 i = 0; i < rec_.titleLen; ++i)
			stats.symbolFreqs[rec_.title[i]]++;

		return;
	}

	// update fields
	//
	uint32 c_field = 0;
	uint32 start_pos = 0;
	uint32 n_field = stats.fields.size();

	uint32 k;
	for (k = 0; k <= rec_.titleLen && c_field < n_field; ++k)
	{
		stats.symbolFreqs[rec_.title[k]] += (k != rec_.titleLen);

		if (rec_.title[k] != stats.fields[c_field].sep && k < rec_.titleLen)
			continue;

		Field& cur_field = stats.fields[c_field];

		// calculate lengths
		//
		if (k - start_pos > cur_field.max_len)
		{
			cur_field.max_len = k - start_pos;
		}
		else if (k - start_pos < cur_field.min_len)
		{
			cur_field.min_len = k - start_pos;
		}

		cur_field.chars.resize(cur_field.max_len);
		uint32 chars_len = MIN(Field::MAX_FIELD_STAT_LEN, k-start_pos);

		// count freqs
		//
		for (uint32 x = 0; x < chars_len; ++x)
		{
			cur_field.chars[x][rec_.title[start_pos+x]]++;
		}
		for (uint32 x = Field::MAX_FIELD_STAT_LEN; x < k-start_pos; ++x)
		{
			cur_field.chars[Field::MAX_FIELD_STAT_LEN][rec_.title[start_pos+x]]++;
		}

		if (cur_field.is_constant)
		{
			if (k - start_pos != cur_field.len)
			{
				cur_field.is_constant = false;
			}
			else
			{
				cur_field.is_constant = std::equal(cur_field.data, cur_field.data + cur_field.len, rec_.title+start_pos);
			}
		}

		if (cur_field.is_len_constant)
		{
			cur_field.is_len_constant = cur_field.len == k-start_pos;
		}

		if (cur_field.is_numeric)
		{
			uint32 value;
			cur_field.is_numeric = core::is_num(rec_.title+start_pos, k-start_pos, value);


			if (cur_field.is_numeric)
			{
				UpdateNumericField(cur_field, value, prevFieldValues[c_field]);

				prevFieldValues[c_field] = value;
			}
		}

		if (!cur_field.is_constant)
		{
			for (uint32 p = 0; p < k - start_pos && p < cur_field.len; ++p)
			{
				cur_field.ham_mask[p] &= cur_field.data[p] == rec_.title[p+start_pos];
			}
		}

		start_pos = k+1;
		c_field++;
	}

	// invalid field
	//
	if ((c_field != n_field) || (k != rec_.titleLen + 1U))
	{
		stats.mixedFormatting = true;
	}

	recordCounter++;
}

void TagAnalyzer::UpdateNumericField(Field &field_, int32 curValue_, int32 prevValue_)
{
	// gather value info
	//
	if (curValue_ < field_.min_value)
	{
		field_.min_value = curValue_;
	}
	else if (curValue_ > field_.max_value)
	{
		field_.max_value = curValue_;
	}

	// gather value RLE stats
	//
	if (recordCounter > 0)
	{
		// get RLE
		//! TODO: checks for ratio
		if (field_.rle_val.cur_sym != curValue_)
		{
			field_.rle_val.run_len++;	// == field_.rle_delta.lens.size()
			field_.rle_val.cur_sym = curValue_;
			field_.rle_val.lens.push_back(field_.rle_val.cur_len);
			field_.rle_val.cur_len = 0;
		}
		else
		{
			field_.rle_val.cur_len++;

			if (field_.rle_val.cur_len > 255)
			{
				field_.rle_val.lens.push_back(255);
				field_.rle_val.cur_len = 0;
				field_.rle_val.run_len++;
			}
		}

		//! TODO: better logic
		//
		if (field_.num_values.size())
		{
			field_.num_values[curValue_]++;
			if (field_.num_values.size() > Field::MAX_NUM_VAL_HUF)
			{
				field_.num_values.clear();
			}
		}
	}
	else
	{
		// clear RLE
		field_.rle_val.cur_sym = curValue_;
		field_.rle_val.cur_len = 0;
		field_.rle_val.run_len = 0;
		field_.rle_val.lens.clear();

		field_.num_values[curValue_]++;
	}


	// gather DETLA info
	//
	if (recordCounter >= 1)
	{
		int32 dvalue = curValue_ - prevValue_;
		if (recordCounter > 1)
		{
			if (dvalue > field_.max_delta)
			{
				field_.max_delta = dvalue;
			}
			else if (dvalue < field_.min_delta)
			{
				field_.min_delta = dvalue;
			}

			// get RLE
			//! TODO: checks for ratio
			if (field_.rle_delta.cur_sym != dvalue)
			{
				field_.rle_delta.run_len++;	// == field_.rle_delta.lens.size()
				field_.rle_delta.cur_sym = dvalue;

				field_.rle_delta.lens.push_back(field_.rle_delta.cur_len);
				field_.rle_delta.cur_len = 0;
			}
			else
			{
				field_.rle_delta.cur_len++;

				if (field_.rle_delta.cur_len > 255)
				{
					field_.rle_delta.lens.push_back(255);
					field_.rle_delta.cur_len = 0;
					field_.rle_delta.run_len++;
				}
			}

			//! TODO: better logic
			//
			if (field_.delta_values.size())
			{
				field_.delta_values[dvalue]++;
				if (field_.delta_values.size() > Field::MAX_NUM_VAL_HUF)
				{
					field_.delta_values.clear();
				}
			}
		}
		else // j == 1
		{
			field_.max_delta = dvalue;
			field_.min_delta = dvalue;

			// clear RLE
			field_.rle_delta.cur_sym = dvalue;
			field_.rle_delta.cur_len = 0;
			field_.rle_delta.run_len = 0;
			field_.rle_delta.lens.clear();

			field_.delta_values[dvalue]++;
		}
	}
}

void TagAnalyzer::FinalizeFieldsStats()
{
	if (stats.mixedFormatting)
		return;

	// Find better encoding of numeric values
	//
	for (std::vector<Field>::iterator f = stats.fields.begin(); f != stats.fields.end(); ++f)
	{
		if (!f->is_numeric)
		{
			if (!f->is_constant)
			{
				f->chars.resize(MIN(f->max_len, Field::MAX_FIELD_STAT_LEN+1));
				f->no_of_bits_per_len = core::bit_length(f->max_len - f->min_len);
			}
			continue;
		}

		int32 diff;
		if (f->max_value - f->min_value < f->max_delta - f->min_delta)
		{
			f->is_delta_coding = false;
			diff = f->max_value - f->min_value;
		}
		else
		{
			f->is_delta_coding = true;
			diff = f->max_delta - f->min_delta;
		}

		f->rle_val.lens.push_back(f->rle_val.cur_len);
		if (f->rle_val.cur_len > 0)
		{
			f->rle_val.cur_len = 0;
			f->rle_val.run_len++;
		}

		if ( (float)recordCounter / (f->rle_val.run_len) > 1.25f)
		{
			f->try_rle_val = true;
		}


		if (f->is_delta_coding)
		{
			f->is_delta_const = diff == 0;

			// decide about RLE
			if (!f->is_delta_const)
			{
				f->rle_delta.lens.push_back(f->rle_delta.cur_len);
				if (f->rle_delta.cur_len > 0)
				{
					f->rle_delta.cur_len = 0;
					f->rle_delta.run_len++;
				}

				if ( (float)recordCounter / (f->rle_delta.run_len) > 1.25f)
				{
					f->try_rle_delta = true;
				}
			}
		}

		// select scheme
		// still remains open a comparison between ValueRle and DeltaRle
		if (f->is_delta_coding && f->is_delta_const)
			f->numeric_scheme = Field::DeltaConst;
		else if (f->is_delta_coding && f->try_rle_delta)
			f->numeric_scheme = Field::DeltaRle;
		else if (f->try_rle_val)
			f->numeric_scheme = Field::ValueRle;
		else if (f->is_delta_coding)
		{
			f->numeric_scheme = Field::DeltaVar;
			uint32 diff = (uint32)(f->max_delta - f->min_delta) + 1;
			f->var_stat_encode = diff <= (int32)Field::MAX_NUM_VAL_HUF && f->delta_values.size();
		}
		else
		{
			f->numeric_scheme = Field::ValueVar;
			uint32 diff = (uint32)(f->max_value - f->min_value) + 1;
			f->var_stat_encode = diff <= (int32)Field::MAX_NUM_VAL_HUF && f->num_values.size();
		}

		f->no_of_bits_per_num = core::bit_length(diff);
		diff = f->max_value - f->min_value;
		f->no_of_bits_per_value = core::bit_length(diff);
	}
}




void TagTokenizerEncoder::StartEncoding(BitMemoryWriter& writer_, TagStats* stats_)
{
	ASSERT(stats == NULL);
	ASSERT(stats_ != NULL);
	ASSERT(!stats_->mixedFormatting);
	stats = stats_;

	StoreFields(writer_);

	recordCounter = 0;
	prevFieldValues.resize(stats->fields.size(), 0);
}

void TagTokenizerEncoder::StoreFields(BitMemoryWriter &bit_stream)
{
	ASSERT(!stats->mixedFormatting);

	byte fieldCount = stats->fields.size();
	bit_stream.PutByte(fieldCount);

	for (std::vector<Field>::iterator f = stats->fields.begin(); f != stats->fields.end(); ++f)
	{
		bit_stream.PutByte(f->sep);				// <---+ put on flags (low 6 bits)
		bit_stream.PutByte(f->is_constant);		// <-- use flags
		if (f->is_constant)
		{
			bit_stream.PutWord(f->len);			// <-- can use 16 bits
			bit_stream.PutBytes(f->data, f->len);
			continue;
		}

		bit_stream.PutByte(f->is_numeric);			// <-- use flags
		if (f->is_numeric)
		{
			bit_stream.PutByte(f->numeric_scheme);	// <-- use flags

			bit_stream.PutWord(f->min_value);		// requred for computing bits_per_value
			bit_stream.PutWord(f->max_value);

			switch (f->numeric_scheme)
			{
			case Field::DeltaConst:
			case Field::DeltaRle:
			case Field::DeltaVar:
				bit_stream.PutWord(f->min_delta);
				bit_stream.PutWord(f->max_delta);

				if (f->numeric_scheme == Field::DeltaVar)
				{
					bit_stream.PutByte((byte)f->var_stat_encode);	// <-- use flags

					if (f->var_stat_encode)							// few values, so use Huffman for them
					{
						uint32 diff = (uint32)(f->max_delta - f->min_delta);
						diff++;

						HuffmanEncoder* huf = f->Huffman_global = new HuffmanEncoder(Field::HUF_GLOBAL_SIZE);
						for (uint32 j = 0; j < diff; ++j)
						{
							huf->Insert(f->delta_values[f->min_delta + j]);
						}
						huf->Complete();
						huf->StoreTree(bit_stream);
					}
				}
				break;

			case Field::ValueRle:
				break;

			case Field::ValueVar:
				{
					bit_stream.PutByte((byte)f->var_stat_encode);	// <-- use flags
					if (f->var_stat_encode)							// few values, so use Huffman for them
					{
						uint32 diff = (uint32)(f->max_value - f->min_value);
						diff++;

						HuffmanEncoder* huf = f->Huffman_global = new HuffmanEncoder(Field::HUF_GLOBAL_SIZE);
						for (uint32 j = 0; j < diff; ++j)
						{
							huf->Insert(f->num_values[f->min_value + j]);
						}
						huf->Complete();
						huf->StoreTree(bit_stream);

					}
					break;
				}

			case Field::None:
				ASSERT(0);
				break;
			}

			continue;
		}

		bit_stream.PutByte(f->is_len_constant);			// <-- use flags
		bit_stream.PutWord(f->len);						// <--
		bit_stream.PutWord(f->max_len);					// <--+ use 16 bits
		bit_stream.PutWord(f->min_len);					// <--
		bit_stream.PutBytes(f->data, f->len);

		for (uint32 j = 0; j < f->len; ++j)
		{
			bit_stream.PutBit(f->ham_mask[j]);
		}
		bit_stream.FlushPartialWordBuffer();

		f->Huffman_local.resize(MIN(f->max_len+1, Field::MAX_FIELD_STAT_LEN+1));

		for (uint32 j = 0; j < MIN(f->max_len, Field::MAX_FIELD_STAT_LEN); ++j)
		{
			f->Huffman_local[j] = NULL;
			if (j >= f->len || !f->ham_mask[j])
			{
				HuffmanEncoder* huf = f->Huffman_local[j] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
				for (uint32 k = 0; k < Field::HUF_LOCAL_SIZE; ++k)
				{
					huf->Insert(f->chars[j][k]);
				}
				huf->Complete(true);
				huf->StoreTree(bit_stream);
			}
		}
		if (f->max_len >= Field::MAX_FIELD_STAT_LEN)
		{
			HuffmanEncoder* huf = f->Huffman_local[Field::MAX_FIELD_STAT_LEN] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
			for (uint32 k = 0; k < Field::HUF_LOCAL_SIZE; ++k)
			{
				huf->Insert(f->chars[(uchar) Field::MAX_FIELD_STAT_LEN][k]);
			}
			huf->Complete(true);
			huf->StoreTree(bit_stream);
		}
	}
}

void TagTokenizerEncoder::EncodeNextFields(BitMemoryWriter &bit_stream, const FastqRecord &rec_)
{
	ASSERT(stats != NULL);

	uint32 c_field = 0;
	uint32 start_pos = 0;

	// store title
	//
	for (uint32 k = 0; k <= rec_.titleLen; ++k)
	{
		ASSERT(c_field < stats->fields.size());
		Field &cur_field = stats->fields[c_field];

		if (rec_.title[k] != cur_field.sep && k < rec_.titleLen)
			continue;

		if (cur_field.is_constant)
		{
			start_pos = k+1;
			c_field++;
			continue;
		}

		if (cur_field.is_numeric)
		{
			int32 value = core::to_num(rec_.title + start_pos, k-start_pos);

			StoreNumericField(bit_stream, cur_field, value, prevFieldValues[c_field]);

			prevFieldValues[c_field] = value;
			start_pos = k+1;
			c_field++;
			continue;
		}

		if (!cur_field.is_len_constant)
		{
			bit_stream.PutBits(k-start_pos - cur_field.min_len, cur_field.no_of_bits_per_len);
		}

		for (uint32 j = 0; j < k-start_pos; ++j)
		{
			if (j >= cur_field.len || !cur_field.ham_mask[j])
			{
				uchar c = rec_.title[start_pos+j];
				const HuffmanEncoder::Code* codes = cur_field.Huffman_local[MIN(j, Field::MAX_FIELD_STAT_LEN)]->GetCodes();
				bit_stream.PutBits(codes[c].code, codes[c].len);
			}
		}

		start_pos = k+1;
		c_field++;
	}

	recordCounter++;
}

void TagTokenizerEncoder::StoreNumericField(BitMemoryWriter &bit_stream, Field &field_, int32 curValue_, int32 prevValue_)
{
	if (recordCounter == 0)
	{
		int32 dval = curValue_ - field_.min_value;
		bit_stream.PutBits(dval, field_.no_of_bits_per_value);

		if (field_.numeric_scheme == Field::ValueRle)
		{
			field_.rle_val.run_len = 0;
			field_.rle_val.cur_len = field_.rle_val.lens[field_.rle_val.run_len];
			field_.rle_val.cur_sym = dval;
			bit_stream.PutBits(field_.rle_val.cur_len, 8);
		}
		return;
	}

	switch (field_.numeric_scheme)
	{
	case Field::DeltaConst:
		break;

	case Field::DeltaRle:
		{
			int32 dval = curValue_ - prevValue_ - field_.min_delta;

			if (recordCounter == 1)
			{
				field_.rle_delta.run_len = 0;
				field_.rle_delta.cur_len = field_.rle_delta.lens[field_.rle_delta.run_len];
				field_.rle_delta.cur_sym = dval;

				bit_stream.PutBits(dval, field_.no_of_bits_per_num);
				bit_stream.PutBits(field_.rle_delta.cur_len, 8);			// optimize
			}
			else
			{
				if (field_.rle_delta.cur_len == 0)
				{
					ASSERT(field_.rle_delta.cur_sym != dval
						   || (field_.rle_delta.lens[field_.rle_delta.run_len] == 255 && field_.rle_delta.cur_sym == dval) );

					field_.rle_delta.run_len++;
					field_.rle_delta.cur_len = field_.rle_delta.lens[field_.rle_delta.run_len];
					field_.rle_delta.cur_sym = dval;

					bit_stream.PutBits(dval, field_.no_of_bits_per_num);
					bit_stream.PutBits(field_.rle_delta.cur_len, 8);		// optimize
				}
				else
				{
					ASSERT(field_.rle_delta.cur_sym == dval);
					field_.rle_delta.cur_len--;
				}
			}

			break;
		}

	case Field::DeltaVar:
		{
			int32 to_store = curValue_ - prevValue_ - field_.min_delta;

			if (field_.Huffman_global)
			{
				const HuffmanEncoder::Code* codes = field_.Huffman_global->GetCodes();
				bit_stream.PutBits(codes[to_store].code, codes[to_store].len);
			}
			else
			{
				bit_stream.PutBits(to_store, field_.no_of_bits_per_num);
			}

			break;
		}

	case Field::ValueRle:
		{
			int32 dval = curValue_ - field_.min_value;

			if (field_.rle_val.cur_len == 0)
			{
				ASSERT(field_.rle_val.cur_sym != dval
					   || (field_.rle_val.lens[field_.rle_val.run_len] == 255 && field_.rle_val.cur_sym == dval) );

				field_.rle_val.run_len++;
				field_.rle_val.cur_len = field_.rle_val.lens[field_.rle_val.run_len];
				field_.rle_val.cur_sym = dval;

				bit_stream.PutBits(dval, field_.no_of_bits_per_value);
				bit_stream.PutBits(field_.rle_val.cur_len, 8);		// optimize
			}
			else
			{
				ASSERT(field_.rle_val.cur_sym == dval);
				field_.rle_val.cur_len--;
			}

			break;
		}

	case Field::ValueVar:
		{
			int32 to_store = curValue_ - field_.min_value;

			if (field_.Huffman_global)
			{
				const HuffmanEncoder::Code* codes = field_.Huffman_global->GetCodes();
				bit_stream.PutBits(codes[to_store].code, codes[to_store].len);
			}
			else
			{
				bit_stream.PutBits(to_store, field_.no_of_bits_per_num);
			}
			break;
		}

	case Field::None:
		ASSERT(0);
		break;
	}
}

void TagTokenizerEncoder::FinishEncoding(BitMemoryWriter &writer_)
{
	ASSERT(stats != NULL);
	ASSERT(!stats->mixedFormatting);

	stats = NULL;

	writer_.FlushPartialWordBuffer();
}


void TagTokenizerDecoder::StartDecoding(BitMemoryReader &reader_)
{
	ReadFields(reader_);

	recordCounter = 0;
	prevFieldValues.resize(stats.fields.size(), 0);
}

void TagTokenizerDecoder::ReadFields(BitMemoryReader &bit_stream)
{
	uint32 n_field = bit_stream.GetByte();
	ASSERT(n_field > 0);

	stats.fields.clear();
	stats.fields.resize(n_field);

	for (uint32 i = 0; i < n_field; ++i)
	{
		Field& field = stats.fields[i];
		field.sep = bit_stream.GetByte();
		field.is_constant = bit_stream.GetByte() != 0;

		if (field.is_constant)
		{
			field.len = bit_stream.GetWord();
			ASSERT(field.len < (1 << 10));
			field.data = new uchar[field.len+1];
			bit_stream.GetBytes(field.data, field.len);
			continue;
		}

		field.is_numeric = bit_stream.GetByte() != 0;
		if (field.is_numeric)
		{
			field.numeric_scheme = bit_stream.GetByte();
			ASSERT(field.numeric_scheme >= Field::ValueVar && field.numeric_scheme <= Field::DeltaConst);

			field.min_value = bit_stream.GetWord();
			field.max_value = bit_stream.GetWord();
			ASSERT(field.min_value <= field.max_value);

			field.no_of_bits_per_value = core::bit_length(field.max_value - field.min_value);
			field.no_of_bits_per_num = 0;

			switch (field.numeric_scheme)
			{
			case Field::DeltaConst:
			case Field::DeltaRle:
			case Field::DeltaVar:
				field.min_delta = (int32)bit_stream.GetWord();
				field.max_delta = (int32)bit_stream.GetWord();
				ASSERT(field.min_delta <= field.max_delta);

				field.no_of_bits_per_num = core::bit_length(field.max_delta - field.min_delta);

				field.is_delta_coding = true;
				field.is_delta_const = field.numeric_scheme == Field::DeltaConst;

				if (field.numeric_scheme == Field::DeltaVar)
				{
					field.var_stat_encode = bit_stream.GetByte();
					if (field.var_stat_encode)			// few values, so use Huffman for them
					{
						ASSERT(field.Huffman_global == NULL);
						field.Huffman_global = new HuffmanEncoder();
						field.Huffman_global->LoadTree(bit_stream);
					}
				}
				break;

			case Field::ValueRle:
				field.no_of_bits_per_num = field.no_of_bits_per_value;
				field.is_delta_coding = false;
				break;

			case Field::ValueVar:
				{
					field.no_of_bits_per_num = field.no_of_bits_per_value;

					field.var_stat_encode = bit_stream.GetByte();
					if (field.var_stat_encode)			// few values, so use Huffman for them
					{
						ASSERT(field.Huffman_global == NULL);
						field.Huffman_global = new HuffmanEncoder();
						field.Huffman_global->LoadTree(bit_stream);
					}
					field.is_delta_coding = false;
					break;
				}

			case Field::None:
			default:
				ASSERT(0);
				break;
			}

			continue;
		}

		field.is_len_constant = bit_stream.GetByte() != 0;
		field.len = bit_stream.GetWord();
		ASSERT(field.len < (1 << 10));

		field.max_len = bit_stream.GetWord();
		ASSERT(field.max_len < (1 << 10));

		field.min_len = bit_stream.GetWord();
		ASSERT(field.min_len < (1 << 10));

		field.no_of_bits_per_len = core::bit_length(field.max_len - field.min_len);
		field.data = new uchar[field.len+1];
		bit_stream.GetBytes(field.data, field.len);
		field.ham_mask = new bool[field.len+1];

		for (uint32 j = 0; j < field.len; ++j)
		{
			field.ham_mask[j] = bit_stream.GetBit() != 0;
		}
		bit_stream.FlushInputWordBuffer();

		field.Huffman_local.resize(MIN(field.max_len, Field::MAX_FIELD_STAT_LEN+1));

		for (uint32 j = 0; j < MIN(field.max_len, Field::MAX_FIELD_STAT_LEN); ++j)
		{
			field.Huffman_local[j] = NULL;
			if (j >= field.len || !field.ham_mask[j])
			{
				field.Huffman_local[j] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
				field.Huffman_local[j]->LoadTree(bit_stream);
			}
		}
		if (field.max_len >= Field::MAX_FIELD_STAT_LEN)
		{
			field.Huffman_local[Field::MAX_FIELD_STAT_LEN] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
			field.Huffman_local[Field::MAX_FIELD_STAT_LEN]->LoadTree(bit_stream);
		}
	}
}

void TagTokenizerDecoder::DecodeNextFields(BitMemoryReader &bit_stream, FastqRecord &rec_)
{
	uint32 n_fields = stats.fields.size();

	for (uint32 j = 0; j < n_fields; ++j)
	{
		Field &cur_field = stats.fields[j];
		if (cur_field.is_constant)
		{
			std::copy(cur_field.data, cur_field.data + cur_field.len, rec_.title + rec_.titleLen);
			rec_.titleLen += cur_field.len;
			rec_.title[rec_.titleLen++] = cur_field.sep;
			continue;
		}
		if (cur_field.is_numeric)
		{
			uint32 num_val = ReadNumericField(bit_stream, cur_field, prevFieldValues[j]);

			rec_.titleLen += core::to_string(rec_.title + rec_.titleLen, num_val);
			prevFieldValues[j] = num_val;

			//rec_.AppendTitle(cur_field.sep);
			rec_.title[rec_.titleLen++] = cur_field.sep;

			continue;
		}

		uint32 field_len = 0;
		if (!cur_field.is_len_constant)
		{
			field_len = bit_stream.GetBits(cur_field.no_of_bits_per_len);
			field_len += cur_field.min_len;
		}
		else
		{
			field_len = cur_field.len;
		}
		ASSERT(field_len <= cur_field.max_len);

		for (uint32 k = 0; k < field_len; ++k)
		{
			if (k < cur_field.len && cur_field.ham_mask[k])
			{
				//rec_.AppendTitle(cur_field.data[k]);
				rec_.title[rec_.titleLen++] = cur_field.data[k];
			}
			else
			{
				HuffmanEncoder *cur_huf = cur_field.Huffman_local[MIN(k, Field::MAX_FIELD_STAT_LEN)];

				uint32 bit = bit_stream.GetBits(cur_huf->GetMinLen());
				int32 h_tmp = cur_huf->DecodeFast(bit);

				while (h_tmp < 0)
				{
					bit = bit_stream.GetBit();
					h_tmp = cur_huf->Decode(bit);
				};

				rec_.title[rec_.titleLen++] = h_tmp;
			}
		}

		rec_.title[rec_.titleLen++] = cur_field.sep;
	}

	rec_.titleLen--;			// do not count last separator

	recordCounter++;
}

uint32 TagTokenizerDecoder::ReadNumericField(BitMemoryReader &bit_stream, Field &field_, int32 prevValue_)
{
	uint32 num_val = 0;
	if (recordCounter == 0)
	{
		num_val = bit_stream.GetBits(field_.no_of_bits_per_value);

		if (field_.numeric_scheme == Field::ValueRle)
		{
			field_.rle_delta.cur_len = bit_stream.GetBits(8);
			field_.rle_delta.cur_sym = (int32)num_val;
		}

		num_val += field_.min_value;
	}
	else
	{
		switch (field_.numeric_scheme)
		{
		case Field::DeltaConst:	// num val = 0
			num_val += prevValue_ + field_.min_delta;
			break;

		case Field::DeltaRle:
			{
				if (recordCounter == 1)
				{
					num_val = bit_stream.GetBits(field_.no_of_bits_per_num);
					field_.rle_delta.cur_sym = (int32)num_val;
					field_.rle_delta.cur_len = bit_stream.GetBits(8);
				}
				else
				{
					if (field_.rle_delta.cur_len == 0)
					{
						num_val = bit_stream.GetBits(field_.no_of_bits_per_num);
						field_.rle_delta.cur_sym = (int32)num_val;
						field_.rle_delta.cur_len = bit_stream.GetBits(8);
					}
					else
					{
						field_.rle_delta.cur_len--;
						num_val = field_.rle_delta.cur_sym;
					}
				}

				num_val += prevValue_ + field_.min_delta;

				break;
			}

		case Field::ValueVar:
		case Field::DeltaVar:
			{
				if (field_.Huffman_global)
				{
					uint32 bit = bit_stream.GetBits(field_.Huffman_global->GetMinLen());
					int32 h_tmp = field_.Huffman_global->DecodeFast(bit);

					while (h_tmp < 0)
					{
						bit = bit_stream.GetBit();
						h_tmp = field_.Huffman_global->Decode(bit);
					};

					num_val = h_tmp;
				}
				else
				{
					num_val = bit_stream.GetBits(field_.no_of_bits_per_num);
				}

				if (field_.numeric_scheme == Field::DeltaVar)
				{
					num_val += prevValue_ + field_.min_delta;
				}
				else
				{
					num_val += field_.min_value;
				}

				break;
			}

		case Field::ValueRle:
			{
				if (field_.rle_delta.cur_len == 0)
				{
					num_val = bit_stream.GetBits(field_.no_of_bits_per_num);
					field_.rle_delta.cur_sym = (int32)num_val;
					field_.rle_delta.cur_len = bit_stream.GetBits(8);
				}
				else
				{
					field_.rle_delta.cur_len--;
					num_val = field_.rle_delta.cur_sym;
				}
			}

			num_val += field_.min_value;

			break;


		case Field::None:
			ASSERT(0);
			break;
		}
	}

	return num_val;
}

void TagTokenizerDecoder::FinishDecoding(BitMemoryReader &reader_)
{
	reader_.FlushInputWordBuffer();
}



void TagRawEncoder::StartEncoding(BitMemoryWriter& writer_, TagStats* stats_)
{
	ASSERT(stats == NULL);
	ASSERT(stats_->mixedFormatting);
	ASSERT(encoder == NULL);

	stats = stats_;

	titleLenBits = core::bit_length(stats->maxTitleLen - stats->minTitleLen);
	writer_.PutWord(stats->minTitleLen);
	writer_.PutWord(stats->maxTitleLen);


	// calculate stats and create huffman tree
	//
	std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
	symbolCount = 0;

	encoder = new HuffmanEncoder(MaxSymbolCount);

	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		if (stats->symbolFreqs[i] > 0)
		{
			symbols[i] = symbolCount++;
			encoder->Insert(stats->symbolFreqs[i]);
		}
	}

	encoder->Complete();


	// store symbols
	//
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		writer_.PutBit(symbols[i] != EmptySymbol);
	}
	writer_.FlushPartialWordBuffer();

	encoder->StoreTree(writer_);
}

void TagRawEncoder::EncodeNextFields(BitMemoryWriter &writer_, const FastqRecord &rec_)
{
	ASSERT(encoder != NULL);

	if (titleLenBits > 0)
		writer_.PutBits(rec_.titleLen - stats->minTitleLen, titleLenBits);

	const HuffmanEncoder::Code* symbolCodes = encoder->GetCodes();
	for (uint32 i = 0; i < rec_.titleLen; ++i)
	{
		const HuffmanEncoder::Code& code = symbolCodes[(int32)symbols[rec_.title[i]]];
		writer_.PutBits(code.code, code.len);
	}
}

void TagRawEncoder::FinishEncoding(BitMemoryWriter &writer_)
{
	ASSERT(encoder != NULL);

	writer_.FlushPartialWordBuffer();

	delete encoder;
	encoder = NULL;
	stats = NULL;
}



void TagRawDecoder::StartDecoding(BitMemoryReader &reader_)
{
	ASSERT(encoder == NULL);

	// read lengths
	//
	minTitleLen = reader_.GetWord();
	maxTitleLen = reader_.GetWord();
	titleLenBits = core::bit_length(maxTitleLen - minTitleLen);

	// read symbols and huffman tree
	//
	symbolCount = 0;
	std::fill(symbols, symbols + MaxSymbolCount, +EmptySymbol);
	for (uint32 i = 0; i < MaxSymbolCount; ++i)
	{
		if (reader_.GetBit())
		{
			symbols[symbolCount++] = i;
		}
	}

	ASSERT(symbolCount > 0);

	encoder = new HuffmanEncoder(symbolCount);
	encoder->LoadTree(reader_);
}

void TagRawDecoder::DecodeNextFields(BitMemoryReader &reader_, FastqRecord &rec_)
{
	ASSERT(encoder != NULL);

	if (titleLenBits > 0)
		rec_.titleLen = reader_.GetBits(titleLenBits) + minTitleLen;
	else
		rec_.titleLen = maxTitleLen;

	for (uint32 i = 0; i < rec_.titleLen; ++i)
	{
		uint32 bit = reader_.GetBits(encoder->GetMinLen());
		int32 sidx = encoder->DecodeFast(bit);
		while (sidx < 0)
		{
			bit = reader_.GetBit();
			sidx = encoder->Decode(bit);
		};

		rec_.title[i] = symbols[sidx];
	}
}

void TagRawDecoder::FinishDecoding(BitMemoryReader &reader_)
{
	ASSERT(encoder != NULL);

	reader_.FlushInputWordBuffer();

	delete encoder;
	encoder = NULL;
}




} // namespace comp

} // namespace dsrc
