/*
         Copyright (C) 2000-2017 the YAMBO team
               http://www.yambo-code.org
 
  Authors (see AUTHORS file for details): AM
  
  This file is distributed under the terms of the GNU 
  General Public License. You can redistribute it and/or 
  modify it under the terms of the GNU General Public 
  License as published by the Free Software Foundation; 
  either version 2, or (at your option) any later version.
 
  This program is distributed in the hope that it will 
  be useful, but WITHOUT ANY WARRANTY; without even the 
  implied warranty of MERCHANTABILITY or FITNESS FOR A 
  PARTICULAR PURPOSE.  See the GNU General Public License 
  for more details.
 
  You should have received a copy of the GNU General Public 
  License along with this program; if not, write to the Free 
  Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
  MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
*/
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>

#include "c_defs.h"
int F77_FUNC_(c_fprintf, C_FPRINTF)(char *lfmt, char *msg,char *rfmt, char *sfmt)
{
 if (strcmp(lfmt,"r")==0) fprintf(stderr,"\r");
 if (strcmp(lfmt,"n")==0) fprintf(stderr,"\n");
 if (strcmp(lfmt,"nn")==0) fprintf(stderr,"\n\n");
 fprintf(stderr,sfmt,msg);
 if (strcmp(rfmt,"n")==0) fprintf(stderr,"\n");
 if (strcmp(rfmt,"nn")==0) fprintf(stderr,"\n\n");
 fflush(stderr);
 return 0;
};
