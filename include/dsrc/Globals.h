/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_GLOBALS
#define H_GLOBALS

#ifndef NDEBUG
#	define DEBUG 1
#endif


// Visual Studio warning supression
//
#if defined (_WIN32)
#	define _CRT_SECURE_NO_WARNINGS
#	pragma warning(disable : 4996) // D_SCL_SECURE
#	pragma warning(disable : 4244) // conversion uint64 to uint32
#	pragma warning(disable : 4267)
#	pragma warning(disable : 4800) // conversion byte to bool
#endif


// assertions
//
#if defined(DEBUG) || defined(_DEBUG)
#	include "assert.h"
#	define ASSERT(x) assert(x)
#else
#	define ASSERT(x) (void)(x)
#endif

#include <string>
#include <stdexcept>

namespace dsrc
{

// basic types
//
typedef char				int8;
typedef unsigned char		uchar, byte, uint8;
typedef short int			int16;
typedef unsigned short int	uint16;
typedef int					int32;
typedef unsigned int		uint32;
typedef long long			int64;
typedef unsigned long long	uint64;


// exceptions
//
class DsrcException : public std::exception
{
	std::string message;

public:
	DsrcException(const char* msg_)
		: message(msg_)
	{}

	DsrcException(const std::string& msg_)
		: message(msg_)
	{}

	~DsrcException() throw()
	{}

	const char* what() const throw()				// for std::exception interface
	{
		return message.c_str();
	}
};


// settings
//

// TODO: handle default parameters
//
#define SET_FIELD(x) (1 << (x))

class FieldMask
{
public:
	FieldMask() : mask(0) {}
	FieldMask(const FieldMask& m_) : mask(m_.mask) {}

	FieldMask AddField(uint32 i_) const
	{
		FieldMask m(*this);
		m.mask |= SET_FIELD(i_);
		return m;
	}

	uint64 GetMask() const
	{
		return mask;
	}

private:
	uint64 mask;
};

struct DsrcCompressionSettings
{
	static const uint32 MinDnaCompressionLevel = 0;
	static const uint32 MaxDnaCompressionLevel = 3;
	static const uint32 MinQualityCompressionLevel = 0;
	static const uint32 MaxQualityCompressionLevel = 3;
	static const uint32 MinFastqBufferSizeMB = 1;
	static const uint32 MaxFastqBufferSizeMB = 1024;
	static const uint32 TagPreserveAllMask = 0;

	static const uint32 DefaultDnaCompressionLevel = 0;
	static const uint32 DefaultQualityCompressionLevel = 0;
	static const uint32 DefaultFastqBufferSizeMB = 8;

	static const bool DefaultLossyQualityCompressionMode = false;
	static const bool DefaultCrc32Calculation = false;

	uint32 dnaCompressionLevel;
	uint32 qualityCompressionLevel;
	uint64 tagPreserveMask;
	bool lossyQualityCompression;
	bool calculateCrc32;
	uint64 fastqBufferSizeMb;

	DsrcCompressionSettings()
		:	dnaCompressionLevel(DefaultDnaCompressionLevel)
		,	qualityCompressionLevel()
		,	tagPreserveMask(TagPreserveAllMask)
		,	lossyQualityCompression(DefaultLossyQualityCompressionMode)
		,	calculateCrc32(DefaultCrc32Calculation)
		,	fastqBufferSizeMb(DefaultFastqBufferSizeMB)
	{}

	static DsrcCompressionSettings Default()
	{
		DsrcCompressionSettings config;
		return config;
	}
};

struct FastqDatasetType
{
	static const uint32 AutoQualityOffsetSelect = 0;
	static const uint32 StandardQualityOffset = 33;

	static const bool DefaultIsColorSpace = false;
	static const bool DefaultHasPlusRepetition = false;

	bool colorSpace;
	bool plusRepetition;
	uint32 qualityOffset;

	FastqDatasetType()
		:	colorSpace(DefaultIsColorSpace)
		,	plusRepetition(DefaultHasPlusRepetition)
		,	qualityOffset(AutoQualityOffsetSelect)
	{}

	static FastqDatasetType Deault()
	{
		FastqDatasetType fq;
		return fq;
	}
};

} // namespace dsrc


#endif // H_GLOBALS
