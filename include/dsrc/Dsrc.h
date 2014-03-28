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

typedef DsrcException DsrcException;
typedef wrap::FastqRecord FastqRecord;
typedef wrap::FastqFile FastqFile;

typedef wrap::FieldMask FieldMask;
typedef wrap::DsrcModule DsrcModule;
typedef wrap::DsrcArchive DsrcArchive;

}

}


#endif // H_DSRC
