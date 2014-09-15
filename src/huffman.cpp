/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "huffman.h"

#include <algorithm>

#include "utils.h"

namespace dsrc
{

namespace comp
{

// ********************************************************************************************
HuffmanEncoder::HuffmanEncoder(uint32 initial_size)
	:	size(initial_size)
	,	n_symbols(0)
	,	min_len(1)
	,	root_id(0)
	,	cur_id(0)
	,	tmp_id(0)
	,	bits_per_id(0)
	,	speedup_tree(NULL)
	,	bit_memory_r(NULL)
	,	bit_memory_w(NULL)
{
	if (size)
	{
		tree  = new Node[2*size-1];
		codes = new Code[2*size-1];
		heap  = new Frequency[size];
	}
	else
	{
		tree  = NULL;
		codes = NULL;
		heap  = NULL;
	}
}

// ********************************************************************************************
HuffmanEncoder::~HuffmanEncoder()
{
	if (tree)
		delete[] tree, tree = NULL;
	if (codes)
		delete[] codes, codes = NULL;
	if (heap)
		delete[] heap, heap = NULL;

	if (speedup_tree)
		delete[] speedup_tree, speedup_tree = NULL;

	ASSERT(bit_memory_w == NULL);
	ASSERT(bit_memory_r == NULL);
}

// ********************************************************************************************
HuffmanEncoder::HuffmanEncoder(const HuffmanEncoder &encoder_)
	:	size(encoder_.size)
	,	n_symbols(0)
	,	min_len(1)
	,	root_id(0)
	,	cur_id(0)
	,	tmp_id(0)
	,	bits_per_id(0)
	,	speedup_tree(NULL)
	,	bit_memory_r(NULL)
	,	bit_memory_w(NULL)
{
	if (size)
	{
		tree  = new Node[2*size-1];
		codes = new Code[2*size-1];
		heap  = new Frequency[size];
	}
	else
	{
		tree  = NULL;
		codes = NULL;
		heap  = NULL;
	}
}

// ********************************************************************************************
HuffmanEncoder::Code* HuffmanEncoder::Complete(bool compact)
{
	if (!n_symbols)
		return NULL;

	// Handle the special case - when maximum no. of symbols == 2
	// but the after compacting the actual the size equals 1 or 0 (*)
	if (n_symbols < 2)
		n_symbols = 2;

	// Make heap of symbols
	std::make_heap(heap, heap+n_symbols);
	
	// Prepare leaves of the tree
	for (uint32 i = 0; i < n_symbols; ++i)
	{
		codes[i].code = 0;
		codes[i].len  = 0;
		tree[i].left_child = -1;
		tree[i].right_child = -1;
	}
	for (uint32 i = n_symbols; i < 2*n_symbols-1; ++i)
	{
		codes[i].code = 0;
		codes[i].len  = 0;
	}

	// Build tree
	int32 heap_size = n_symbols;

	// Remove symbols with 0 frequency
	if (compact)
	{
		// (*) special case support
		if (heap_size == 2 && heap[0].frequency == 0)
		{
			heap[0].frequency = 1;
			if (heap[1].frequency == 0)
				heap[1].frequency = 1;
		}
		else
		{
			while (heap_size > 2 && heap[0].frequency == 0)
				std::pop_heap(heap, heap+heap_size--);
		}
	}

	int32 present_symbols = heap_size;

	if (!present_symbols)
		return codes;

	for (int32 i = 0; i < present_symbols-1; ++i)
	{
		Frequency left = heap[0];
		std::pop_heap(heap, heap+heap_size--);
		Frequency right = heap[0];
		std::pop_heap(heap, heap+heap_size--);
		
		heap[heap_size].symbol = n_symbols+i;
		heap[heap_size].frequency = left.frequency + right.frequency;
		std::push_heap(heap, heap+ ++heap_size);

		tree[n_symbols+i].left_child  = left.symbol;
		tree[n_symbols+i].right_child = right.symbol;
	}

	// Compute codes
	for (int32 i = n_symbols+present_symbols-2; i >= (int32)n_symbols; --i)
	{
		codes[tree[i].left_child].len   = codes[i].len+1;
		codes[tree[i].left_child].code  = (codes[i].code << 1);
		codes[tree[i].right_child].len  = codes[i].len+1;
		codes[tree[i].right_child].code = (codes[i].code << 1) | 1;
	}

	root_id = n_symbols + present_symbols - 2;
	cur_id = root_id;

	return codes;
} 

// ********************************************************************************************
void HuffmanEncoder::StoreTree(core::BitMemoryWriter &bit_stream)
{
	ASSERT(n_symbols > 1);

	bit_memory_w = &bit_stream;

	// save data
	//
	bit_memory_w->FlushPartialWordBuffer();

	uint64 sizePos = bit_memory_w->Position();
	bit_memory_w->PutWord(0);				// fill dummy value

	const uint32 tree_size = n_symbols;		// in prev version 'size' was used
	bits_per_id = core::int_log(tree_size, 2);
	if (tree_size & (tree_size-1))			// size is not power of 2
		bits_per_id++;

	ASSERT(bits_per_id > 0);

	min_len = n_symbols;					 // was 32 here and failed for n_symbols == 1
	for (uint32 i = 0; i < n_symbols; ++i)
	{
		if (codes[i].len < min_len && codes[i].len > 0)
			min_len = codes[i].len;
	}

	int32 node_id = root_id;
	bit_memory_w->PutWord(root_id);
	bit_memory_w->PutWord(n_symbols);
	bit_memory_w->PutByte((uchar) min_len);
	EncodeProcess(node_id);
	bit_memory_w->FlushPartialWordBuffer();
	

	// save size
	//
	uint32 memSize = (uint32) bit_memory_w->Position() - sizePos;
	ASSERT(memSize > 0 && memSize < (1 << 20));
	bit_memory_w->SetPosition(sizePos);
	bit_memory_w->PutWord(memSize);
	bit_memory_w->SetPosition(sizePos + memSize);

	bit_memory_w = NULL;
}


// ********************************************************************************************
void HuffmanEncoder::LoadTree(core::BitMemoryReader &bit_stream)
{
	bit_memory_r = &bit_stream;

	bit_memory_r->FlushInputWordBuffer();

	uint32 memBegin = bit_memory_r->Position();
	uint32 memSize = bit_memory_r->GetWord();
	ASSERT(memSize > 0 && memSize < (1 << 20));		// modified

	root_id = bit_memory_r->GetWord();
	n_symbols = bit_memory_r->GetWord();
	ASSERT(n_symbols > 1);
	ASSERT(n_symbols < (1 << 10));

	tmp_id = root_id;
	cur_id = root_id;

	min_len = bit_memory_r->GetByte();
	RestartDecompress(n_symbols, root_id);
	bits_per_id = core::int_log(size, 2);

	if (n_symbols & (n_symbols-1))			// size is not power of 2
		bits_per_id++;

	ASSERT(bits_per_id > 0);

//	n_symbols = tmp;
	int32 node_id = root_id - n_symbols + 1;
	root_id = node_id;
	tmp_id = root_id;
	cur_id = root_id;
	DecodeProcess(node_id);

	bit_memory_r->FlushInputWordBuffer();

	if (!min_len)
		min_len = 1;
	ComputeSpeedupTree();

	ASSERT(memBegin + memSize == bit_memory_r->Position());
	bit_memory_r = NULL;
}

// ********************************************************************************************
void HuffmanEncoder::ComputeSpeedupTree()
{
	if (!min_len)
		return;

	if (speedup_tree)
		delete[] speedup_tree;
	speedup_tree = new int32[(uint32) (1 << min_len)];

	for (int32 i = 0; i < (1 << min_len); ++i)
	{
		cur_id = root_id;
		for (int32 j = min_len-1; j >= 0; --j)
		{
			Decode(i & (1 << j));
		}
		speedup_tree[i] = cur_id;
	}

	cur_id = root_id;
	tmp_id = root_id;
}

// ********************************************************************************************
void HuffmanEncoder::Restart(uint32 _size)
{
	if (size != _size)
	{
		size = _size;

		if (tree)
			delete[] tree;
		if (codes)
			delete[] codes;
		if (heap)
			delete[] heap;
		if (speedup_tree)
			delete[] speedup_tree;

		if (size)
		{
			int* tmp = new int[100];
			codes = new Code[2*size-1];
			heap  = new Frequency[size];
			tree  = new Node[2*size-1];
			delete tmp;
		}
		else
		{
			tree  = NULL;
			codes = NULL;
			heap  = NULL;
		}
		speedup_tree = NULL;
	}

	n_symbols = 0;
}

// ********************************************************************************************
void HuffmanEncoder::RestartDecompress(uint32 _size, uint32 _root_id)
{
	if (tree)
		delete[] tree;
	if (codes)
		delete[] codes;
	if (heap)
		delete[] heap;
	if (speedup_tree)
		delete[] speedup_tree;


	size = _size;
	if (size)
	{
		tree  = new Node[_root_id - _size + 2];
		//		codes = new Code[2*size-1];
		codes = NULL;
		//		heap  = new Frequency[size];
		heap = NULL;
	}
	else
	{
		tree  = NULL;
		codes = NULL;
		heap  = NULL;
	}
	n_symbols = _size;

	speedup_tree = NULL;
}

} // namespace comp

} // namespace dsrc
