/*---------------------------------------------------------------------------*\

  FILE........: gdb_stdio.h
  AUTHOR......: David Rowe
  DATE CREATED: April 23 2013

  Some stdio I/O functions that perform I/O on the host using gdb.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2013 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __GDB_STDIO__
#define __GDB_STDIO__

#include <stdio.h>
#include <stdarg.h>

void gdb_stdio_fprintf(FILE *file, const char *format, ...);
void gdb_stdio_printf(const char *format, ...);
FILE *gdb_stdio_fopen(char file_name[], char mode[]);
void gdb_stdio_fclose(FILE *file);
int gdb_stdio_fwrite(void *ptr, int size, int nmemb, FILE *file);
int gdb_stdio_fread(void *ptr, int size, int nmemb, FILE *file);

#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite

#endif
