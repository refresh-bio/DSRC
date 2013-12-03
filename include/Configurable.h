/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_CONFIGURABLE
#define H_CONFIGURABLE

#include "Globals.h"

namespace dsrc
{

namespace wrap
{

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

class Configurable
{
public:
	Configurable();
	virtual ~Configurable();

	void SetFastqBufferSizeMB(uint64 size_);
	uint64 GetFastqBufferSizeMB() const;

	void SetDnaCompressionLevel(uint32 level_);
	uint32 GetDnaCompressionLevel() const;

	void SetQualityCompressionLevel(uint32 level_);
	uint32 GetQualityCompressionLevel() const;

	void SetLossyCompression(bool lossy_);
	bool IsLossyCompression() const;

	void SetQualityOffset(uint32 off_);
	uint32 GetQualityOffset() const;

	void SetColorSpace(bool cs_);
	bool IsColorSpace() const;

	void SetPlusRepetition(bool rep_);
	bool IsPlusRepetition() const;

	void SetThreadsNumber(uint32 threadNum_);
	uint32 GetThreadsNumber() const;

	void SetStdIoUsing(bool use_);
	bool IsStdIoUsing() const;

	void SetCrc32Checking(bool use_);
	bool IsCrc32Checking() const;

	void SetTagFieldFilterMask(uint64 mask_);
	uint64 GetTagFieldFilterMask() const;

protected:
	void* GetInputParameters() const;

private:
	struct ConfigImpl;
	ConfigImpl* config;
};

} // namespace wrap

} // namespace dsrc

#endif // H_CONFIGURABLE
