/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#include "DsrcFile.h"
#include "DsrcIo.h"
#include "BitMemory.h"

#include <cstring>

namespace dsrc
{

namespace comp
{

using namespace core;

DsrcFileWriter::DsrcFileWriter()
	:	fileStream(NULL)
	,	currentBlockId(0)
{		
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(DsrcFileHeader), 0);
	fileFooter.dummyByte = 0;
}

DsrcFileWriter::~DsrcFileWriter()
{
	if (fileStream != NULL)
		delete fileStream;
}

void DsrcFileWriter::StartCompress(const std::string& fileName_)
{
	ASSERT(fileStream == NULL);

	fileStream = new FileStreamWriterExt(fileName_);

	// clear header and footer
	//
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(DsrcFileHeader), 0);

	fileFooter.blockSizes.clear();
	fileFooter.dummyByte = 0;


	// skip header pos
	//
	fileStream->SetPosition(DsrcFileHeader::HeaderSize);

	currentBlockId = 0;
}

void DsrcFileWriter::WriteNextChunk(const DsrcDataChunk* block_)
{
	ASSERT(block_ != NULL);
	ASSERT(block_->size > 0);

	fileStream->Write(block_->data.Pointer(), block_->size);
	fileFooter.blockSizes.push_back(block_->size);

	for (uint32 i = 0; i < fq::StreamsInfo::StreamCount; ++i)
	{
		fastqStreamInfo.sizes[i] += block_->rawStreamsInfo.sizes[i];
		dsrcStreamInfo.sizes[i] += block_->compStreamsInfo.sizes[i];
	}
	currentBlockId++;
}

void DsrcFileWriter::FinishCompress()
{
	ASSERT(fileFooter.blockSizes.size() > 0);
	ASSERT(fileFooter.blockSizes.size() == currentBlockId);

	// prepare header
	//
	fileHeader.dummyByte = DsrcFileHeader::DummyByteValue;
	fileHeader.versionMajor = DsrcFileHeader::VersionMajor;
	fileHeader.versionMinor = DsrcFileHeader::VersionMinor;
	fileHeader.versionRev = DsrcFileHeader::VersionRev;
	std::fill(fileHeader.reserved, fileHeader.reserved + DsrcFileHeader::ReservedBytes, +DsrcFileHeader::DummyByteValue);
	fileHeader.blockCount = fileFooter.blockSizes.size();
	fileHeader.recordsCount = 0; //! TODO!!!
	fileHeader.footerOffset = fileStream->Position();

	// write footer
	//
	fileFooter.dummyByte = DsrcFileFooter::DummyByteValue;
	WriteFileFooter();

	// fill header and write
	//
	fileHeader.footerSize = fileStream->Position() - fileHeader.footerOffset;

	fileStream->SetPosition(0);
	WriteFileHeader();

	// cool, exit
	//
	fileStream->Close();

	delete fileStream;

	fileStream = NULL;
}

void DsrcFileWriter::WriteFileHeader()
{
	// TODO: here we can just directly flush whole header structure to IO
	//
	BitMemoryWriter writer(DsrcFileHeader::HeaderSize);
	writer.PutByte(fileHeader.dummyByte);
	writer.PutByte(fileHeader.versionMajor);
	writer.PutByte(fileHeader.versionMinor);
	writer.PutByte(fileHeader.versionRev);

	writer.PutWord(fileHeader.footerSize);
	writer.PutDWord(fileHeader.footerOffset);

	writer.PutDWord(fileHeader.recordsCount);
	writer.PutDWord(fileHeader.blockCount);

	writer.PutBytes(fileHeader.reserved, DsrcFileHeader::ReservedBytes);

	fileStream->Write(writer.Pointer(), writer.Position());
}

void DsrcFileWriter::WriteFileFooter()
{	
	// store data
	//
	BitMemoryWriter writer(fileFooter.blockSizes.size() * 4 + DsrcFileFooter::DatasetTypeSize + DsrcFileFooter::CompressionSettingsSize);
	writer.PutByte(fileFooter.dummyByte);

	// store blocks
	//
	writer.PutBytes((byte*)fileFooter.blockSizes.data(), fileFooter.blockSizes.size() * 4);

	// store dataset info
	//
	byte flags = 0;
	if (fileFooter.datasetType.colorSpace)
		flags |= DsrcFileFooter::FLAG_COLOR_SPACE;
	if (fileFooter.datasetType.plusRepetition)
		flags |= DsrcFileFooter::FLAG_PLUS_REPETITION;

	writer.PutByte(flags);
	writer.PutByte(fileFooter.datasetType.qualityOffset);

	// store compression info
	//
	flags = 0;
	if (fileFooter.compSettings.lossy)
		flags |= DsrcFileFooter::FLAG_LOSSY_QUALITY;
	if (fileFooter.compSettings.calculateCrc32)
		flags |= DsrcFileFooter::FLAG_CALCULATE_CRC32;
	writer.PutByte(flags);
	writer.PutByte(fileFooter.compSettings.dnaOrder);
	writer.PutByte(fileFooter.compSettings.qualityOrder);
	writer.PutDWord(fileFooter.compSettings.tagPreserveFlags);

	// flush
	//
	fileStream->Write(writer.Pointer(), writer.Position());
}

DsrcFileReader::DsrcFileReader()
	:	fileStream(NULL)
	,	currentBlockId(0)
{
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(DsrcFileHeader), 0);
	fileFooter.dummyByte = 0;
}

DsrcFileReader::~DsrcFileReader()
{
	if (fileStream)
		delete fileStream;
}

void DsrcFileReader::StartDecompress(const std::string& fileName_)
{
	ASSERT(fileStream == NULL);

	fileStream = new FileStreamReaderExt(fileName_);

	if (fileStream->Size() == 0)
		throw DsrcException("Empty file.");

	// Read file header
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(DsrcFileHeader), 0);
	ReadFileHeader();

	// Check version compatibility
	if (!(fileHeader.versionMajor == DsrcFileHeader::VersionMajor && fileHeader.versionMinor == DsrcFileHeader::VersionMinor))
	{
		//! TODO: add old file version checking
		delete fileStream;
		fileStream = NULL;
		throw DsrcException("Invalid archive or old unsupported version");
	}

	if ((fileHeader.blockCount == 0ULL) || fileHeader.footerOffset + (uint64)fileHeader.footerSize > fileStream->Size())
	{
		delete fileStream;
		fileStream = NULL;
		throw DsrcException("Corrupted DSRC archive header");
	}

	// clean footer
	//
	fileFooter.blockSizes.clear();
	fileFooter.blockSizes.resize(fileHeader.blockCount, 0);

	fileStream->SetPosition(fileHeader.footerOffset);
	ReadFileFooter();

	if (fileFooter.dummyByte != DsrcFileFooter::DummyByteValue)
	{
		delete fileStream;
		fileStream = NULL;
		throw DsrcException("Corrupted DSRC archive footer");
	}

	fileStream->SetPosition(DsrcFileHeader::HeaderSize);

	currentBlockId = 0;
}

bool DsrcFileReader::ReadNextChunk(DsrcDataChunk* block_)
{
	ASSERT(block_ != NULL);

	if (currentBlockId == fileHeader.blockCount)
	{
		block_->size = 0;
		return false;
	}

	block_->size = fileFooter.blockSizes[currentBlockId];
	if (block_->data.Size() < block_->size)
	{
		block_->data.Extend(block_->size);
	}
	fileStream->Read(block_->data.Pointer(), block_->size);
	currentBlockId++;

	return true;
}

void DsrcFileReader::FinishDecompress()
{
	fileStream->Close();

	delete fileStream;
	fileStream = NULL;
}

void DsrcFileReader::ReadFileHeader()
{
	// TODO: here we can just read directly whole header structure from IO
	//
	Buffer buffer(DsrcFileHeader::HeaderSize);
	fileStream->Read(buffer.Pointer(), DsrcFileHeader::HeaderSize);

	BitMemoryReader reader(buffer.Pointer(), DsrcFileHeader::HeaderSize);
	fileHeader.dummyByte = reader.GetByte();
	fileHeader.versionMajor = reader.GetByte();
	fileHeader.versionMinor = reader.GetByte();
	fileHeader.versionRev = reader.GetByte();

	fileHeader.footerSize = reader.GetWord();
	fileHeader.footerOffset = reader.GetDWord();

	fileHeader.recordsCount = reader.GetDWord();
	fileHeader.blockCount = reader.GetDWord();

	reader.GetBytes(fileHeader.reserved, DsrcFileHeader::ReservedBytes);
}

void DsrcFileReader::ReadFileFooter()
{
	Buffer buffer(fileHeader.footerSize);
	fileStream->Read(buffer.Pointer(), fileHeader.footerSize);

	BitMemoryReader reader(buffer.Pointer(), fileHeader.footerSize);
	fileFooter.dummyByte = reader.GetByte();

	// read blocks
	//
	reader.GetBytes((byte*)fileFooter.blockSizes.data(), fileFooter.blockSizes.size() * 4);


	// read dataset info
	//
	byte flags = reader.GetByte();
	fileFooter.datasetType.colorSpace = (flags & DsrcFileFooter::FLAG_COLOR_SPACE) != 0;
	fileFooter.datasetType.plusRepetition = (flags & DsrcFileFooter::FLAG_PLUS_REPETITION) != 0;
	fileFooter.datasetType.qualityOffset = reader.GetByte();

	// read compression info
	//
	flags = reader.GetByte();
	fileFooter.compSettings.lossy = (flags & DsrcFileFooter::FLAG_LOSSY_QUALITY);
	fileFooter.compSettings.calculateCrc32 = (flags & DsrcFileFooter::FLAG_CALCULATE_CRC32);
	fileFooter.compSettings.dnaOrder = reader.GetByte();
	fileFooter.compSettings.qualityOrder = reader.GetByte();
	fileFooter.compSettings.tagPreserveFlags = reader.GetDWord();
}

} // namespace comp

} // namespace dsrc
