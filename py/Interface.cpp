/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include <boost/python/detail/wrap_python.hpp>		// must be included before <Python>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <Python.h>

#include "../include/dsrc/Globals.h"
#include "../include/dsrc/FastqRecord.h"
#include "../include/dsrc/FastqFile.h"
#include "../include/dsrc/DsrcModule.h"
#include "../include/dsrc/DsrcArchive.h"

namespace dsrc
{

namespace py
{

namespace boo = boost::python;
using namespace wrap;

class PyException : public DsrcException
{
public:
	PyException(const char* msg_)
		: DsrcException(msg_)
	{}

	PyException(const std::string& msg_)
		: DsrcException(msg_)
	{}

	static void translate(const PyException& e_)		// boost::python requirement
	{
		PyErr_SetString(PyExc_RuntimeError, e_.what());
	}
};


BOOST_PYTHON_MODULE(pydsrc)
{
	boo::register_exception_translator<PyException>(&PyException::translate);

	boo::class_<FastqRecord>("FastqRecord")
		.def_readwrite("tag", &FastqRecord::tag)
		.def_readwrite("sequence", &FastqRecord::sequence)
		.def_readwrite("plus", &FastqRecord::plus)
		.def_readwrite("quality", &FastqRecord::quality)
	;

	boo::class_<FastqFile, boost::noncopyable>("FastqFile")
		.def("Open", &FastqFile::Open)
		.def("Create", &FastqFile::Create)
		.def("Close", &FastqFile::Close)
		.def("ReadNextRecord", &FastqFile::ReadNextRecord)
		.def("WriteNextRecord", &FastqFile::WriteNextRecord)
	;

	boo::class_<FieldMask>("FieldMask")
		.def("AddField", &FieldMask::AddField)
		.def("GetMask", &FieldMask::GetMask)
	;

	boo::class_<DsrcArchive, boost::noncopyable>("DsrcArchive")
		.def("StartCompress", &DsrcArchive::StartCompress)
		.def("WriteNextRecord", &DsrcArchive::WriteNextRecord)
		.def("FinishCompress", &DsrcArchive::FinishCompress)
		.def("StartDecompress", &DsrcArchive::StartDecompress)
		.def("ReadNextRecord", &DsrcArchive::ReadNextRecord)
		.def("FinishDecompress", &DsrcArchive::FinishDecompress)
		.add_property("LossyCompression", &DsrcArchive::IsLossyCompression, &DsrcArchive::SetLossyCompression)
		.add_property("DNACompressionLevel", &DsrcArchive::GetDnaCompressionLevel, &DsrcArchive::SetDnaCompressionLevel)
		.add_property("QualityCompressionLevel", &DsrcArchive::GetQualityCompressionLevel, &DsrcArchive::SetDnaCompressionLevel)
		.add_property("TagFieldFilterMask", &DsrcArchive::GetTagFieldFilterMask, &DsrcArchive::SetTagFieldFilterMask)
		.add_property("PlusRepetition", &DsrcArchive::IsPlusRepetition, &DsrcArchive::SetPlusRepetition)
		.add_property("QualityOffset", &DsrcArchive::GetQualityOffset, &DsrcArchive::SetQualityOffset)
		.add_property("ColorSpace", &DsrcArchive::IsColorSpace, &DsrcArchive::SetColorSpace)
		.add_property("FastqBufferSizeMB", &DsrcArchive::GetFastqBufferSizeMB, &DsrcArchive::SetFastqBufferSizeMB)
		.add_property("Crc32Checking", &DsrcArchive::IsCrc32Checking, &DsrcArchive::SetCrc32Checking)
	;

	boo::class_<DsrcModule, boost::noncopyable>("DsrcModule")
		.def("Compress", &DsrcModule::Compress)
		.def("Decompress", &DsrcModule::Decompress)
		//.def("Usage", &DsrcModule::Usage)
		.add_property("LossyCompression", &DsrcModule::IsLossyCompression, &DsrcModule::SetLossyCompression)
		.add_property("DNACompressionLevel", &DsrcModule::GetDnaCompressionLevel, &DsrcModule::SetDnaCompressionLevel)
		.add_property("QualityCompressionLevel", &DsrcModule::GetQualityCompressionLevel, &DsrcModule::SetDnaCompressionLevel)
		.add_property("TagFieldFilterMask", &DsrcModule::GetTagFieldFilterMask, &DsrcModule::SetTagFieldFilterMask)
		.add_property("FastqBufferSizeMB", &DsrcModule::GetFastqBufferSizeMB, &DsrcModule::SetFastqBufferSizeMB)
		.add_property("ThreadsNumber", &DsrcModule::GetThreadsNumber, &DsrcModule::SetThreadsNumber)
		.add_property("Crc32Checking", &DsrcModule::IsCrc32Checking, &DsrcModule::SetCrc32Checking)
	;
}

} // namespace py

} // namespace dsrc
