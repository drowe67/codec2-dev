/*---------------------------------------------------------------------------*\

  FILE........: gdb_stdio.c
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

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "gdb_stdio.h"

#define MAX_STR 2048

/* command codes we use to signal host */

#define GDB_STDIO_PRINTF  1
#define GDB_STDIO_FOPEN   2
#define GDB_STDIO_FCLOSE  3
#define GDB_STDIO_FWRITE  4
#define GDB_STDIO_FREAD   5
#define GDB_STDIO_FPRINTF 6

/* globals we use to communicate with host */

volatile int   gdb_stdio_func = 0;
volatile int   gdb_stdio_ret = 0;
volatile char *gdb_stdio_pstr1;
volatile char *gdb_stdio_pstr2;
volatile int   gdb_stdio_strlen1;
volatile int   gdb_stdio_strlen2;
volatile FILE *gdb_stdio_file;
volatile void *gdb_stdio_ptr;
volatile int   gdb_stdio_size;
volatile int   gdb_stdio_nmem;

void gdb_stdio_fprintf(FILE *file, const char *format, ...) {
    va_list arg;
    char str[MAX_STR];

    va_start(arg, format);
    vsnprintf(str, MAX_STR, format, arg);
    va_end(arg);
    gdb_stdio_file = file;
    gdb_stdio_pstr1 = str;
    gdb_stdio_strlen1 = strlen(str);

    gdb_stdio_func = GDB_STDIO_FPRINTF;
    while(gdb_stdio_func);
}

void gdb_stdio_printf(const char *format, ...) {
    va_list arg;
    char str[MAX_STR];

    va_start(arg, format);
    vsnprintf(str, MAX_STR, format, arg);
    va_end(arg);
    gdb_stdio_pstr1 = str;
    gdb_stdio_strlen1 = strlen(str);

    gdb_stdio_func = GDB_STDIO_PRINTF;
    while(gdb_stdio_func);
}

FILE *gdb_stdio_fopen(char file_name[], char mode[]) {
    gdb_stdio_pstr1 = file_name;
    gdb_stdio_pstr2 = mode;
    gdb_stdio_strlen1 = strlen(file_name);
    gdb_stdio_strlen2 = strlen(mode);

    gdb_stdio_func = GDB_STDIO_FOPEN;
    while(gdb_stdio_func);
    return (FILE*)gdb_stdio_ret;
}

void gdb_stdio_fclose(FILE *file) {
    gdb_stdio_file = file;

    gdb_stdio_func = GDB_STDIO_FCLOSE;
    while(gdb_stdio_func);
}

int gdb_stdio_fwrite(void *ptr, int size, int nmem, FILE *file) {
    gdb_stdio_ptr = ptr;
    gdb_stdio_size = size;
    gdb_stdio_nmem = nmem;
    gdb_stdio_file = file;
 
    gdb_stdio_func = GDB_STDIO_FWRITE;
    while(gdb_stdio_func);
    return gdb_stdio_ret;       
}

int gdb_stdio_fread(void *ptr, int size, int nmem, FILE *file) {
    gdb_stdio_ptr = ptr;
    gdb_stdio_size = size;
    gdb_stdio_nmem = nmem;
    gdb_stdio_file = file;
 
    gdb_stdio_func = GDB_STDIO_FREAD;
    while(gdb_stdio_func);
    return gdb_stdio_ret;       
}

