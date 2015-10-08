/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_DSRC
#define H_DSRC

#include "Globals.h"
#include "DsrcModule.h"
#include "FastqFile.h"
#include "DsrcArchive.h"
#include "FastqFile.h"

namespace dsrc
{

namespace lib
{

typedef ext::FastqRecord FastqRecord;
typedef ext::FastqFileRecordsReader FastqFileRecordsReader;
typedef ext::FastqFileRecordsWriter FastqFileRecordsWriter;
typedef ext::FastqFileBlocksReader FastqFileBlocksReader;
typedef ext::FastqFileBlocksWriter FastqFileBlocksWriter;

typedef FieldMask FieldMask;
typedef DsrcCompressionSettings DsrcCompressionSettings;
typedef FastqDatasetType FastqDatasetType;

typedef ext::DsrcModule DsrcModule;


typedef ext::DsrcArchiveRecordsWriter DsrcArchiveRecordsWriter;
typedef ext::DsrcArchiveRecordsReader DsrcArchiveRecordsReader;
typedef ext::DsrcArchiveBlocksWriterST DsrcArchiveBlocksWriter;
typedef ext::DsrcArchiveBlocksReaderST DsrcArchiveBlocksReader;


}



}


#endif // H_DSRC
