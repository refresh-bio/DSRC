/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/


#ifndef _HUFFMAN_H
#define _HUFFMAN_H

#include "../include/dsrc/Globals.h"

#include "BitMemory.h"

namespace dsrc
{

namespace comp
{

// this class is taken almost purely from DSRC 1.0 and needs a-bit-o'-refactoring

class HuffmanEncoder 
{
public:
	struct Code
	{
		uint32 code;
		uint32 len;
	};

	HuffmanEncoder(uint32 initial_size = 0);
	~HuffmanEncoder();
	HuffmanEncoder(const HuffmanEncoder& encoder_);

	void StoreTree(core::BitMemoryWriter& bit_stream);
	void LoadTree(core::BitMemoryReader& bit_stream);

	Code* Complete(bool compact = true);
	inline int32 Decode(const uint32 bit);
	void Restart(uint32 _size = 0);
	void RestartDecompress(uint32 _size = 0, uint32 _root_id = 0);

	inline bool Insert(const uint32 frequency);
	inline int32 DecodeFast(const uint32 bits);

	uint32 GetMinLen() const { return min_len; }
	uint32 GetBitsPerId() const { return bits_per_id; }
	uint32 GerSymbolsNum() const { return n_symbols; }

	const Code* GetCodes() const {return codes;}

private:
	struct Frequency 
	{
		uint32 symbol;
		uint32 frequency;

		Frequency()
			:	symbol(0)
			,	frequency(0)
		{}

		bool operator<(const Frequency &x) const
		{
			return frequency > x.frequency || (frequency == x.frequency && symbol > x.symbol);
		}
	};

	struct Node
	{
		int32 left_child;
		int32 right_child;

		Node()
			:	left_child(0)
			,	right_child(0)
		{}
	};

	uint32 size;
	uint32 n_symbols;
	uint32 min_len;

	int32 root_id;
	int32 cur_id;
	int32 tmp_id;

	uint32 bits_per_id;

	Node *tree;
	Frequency *heap;
	Code *codes;
	int32 *speedup_tree;

	core::BitMemoryReader *bit_memory_r;
	core::BitMemoryWriter *bit_memory_w;

	inline void EncodeProcess(int32 node_id);
	inline int32 DecodeProcess(int32 node_id);

	void ComputeSpeedupTree();
};

// ********************************************************************************************
void HuffmanEncoder::EncodeProcess(int32 node_id)
{
	if (tree[node_id].left_child == -1)		// leaf
	{
		bit_memory_w->PutBit(1);
		bit_memory_w->PutBits(node_id, bits_per_id);
	}
	else
	{
		bit_memory_w->PutBit(0);
		EncodeProcess(tree[node_id].left_child);
		EncodeProcess(tree[node_id].right_child);
	}
}

// ********************************************************************************************
int32 HuffmanEncoder::DecodeProcess(int32 node_id)
{
	ASSERT(node_id >= 0);

	//uint32 flag = bit_memory_r->GetBit();

	//if (!flag)
	if (!bit_memory_r->GetBit())
	{
		--tmp_id;
		tree[node_id].left_child  = DecodeProcess(tmp_id);
		tree[node_id].right_child = DecodeProcess(tmp_id);
		return node_id;
	}
	else
	{
		//uint32 tmp = bit_memory_r->GetBits(bits_per_id);
		//		tree[tmp].left_child  = -1;
		//		tree[tmp].right_child = -1;
		//return -(int32)tmp;
		return -(int32)bit_memory_r->GetBits(bits_per_id);
	}
}

// ********************************************************************************************
bool HuffmanEncoder::Insert(const uint32 frequency)
{
	if (n_symbols == size)
		return false;

	heap[n_symbols].symbol    = n_symbols;
	heap[n_symbols].frequency = frequency;
	n_symbols++;

	return true;
}

// ********************************************************************************************
int32 HuffmanEncoder::Decode(const uint32 bit) 
{
	if (cur_id <= 0)
		cur_id = root_id;
	if (bit)
		cur_id = tree[cur_id].right_child;
	else
		cur_id = tree[cur_id].left_child;

	if (cur_id <= 0)
		return -cur_id;				// Symbol found
	
	return -1;		
}

// ********************************************************************************************
inline int32 HuffmanEncoder::DecodeFast(const uint32 bits)
{
	cur_id = speedup_tree[bits];
//	if (cur_id < n_symbols)
	if (cur_id <= 0)
		return -cur_id;				// Symbol found
	
	return -1;					// Not found yet
}


} // namespace comp

} // namespace dsrc

#endif

