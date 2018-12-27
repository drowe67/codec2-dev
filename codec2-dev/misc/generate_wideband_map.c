/*---------------------------------------------------------------------------*\

  FILE........: generate_wideband_map.c
  AUTHOR......: Phil Ayres
  DATE CREATED: 17 Jul 2017

  Generate header file containing wideband DCT2 map, runs at compile time.
  Adapted from generate_codebook.c

\*---------------------------------------------------------------------------*/

/*
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "c2wideband.h"

static const int Nt = C2WB_NT;
static const int K = C2WB_K;

static const char usage[] =
"Usage: %s filename array_name\n"
"\tCreate C code for wideband DCT2 map.\n";

static const char format[] =
"The table format must be:\n"
"\t8 rows by 30 (Nt x K) floating point numbers to fill the specified dimensions.\n";

static const char header[] =
"/* THIS IS A GENERATED FILE. Edit generate_wideband_map.c and its input */\n\n"
"/*\n"
" * This intermediary file and the files that used to create it are under \n"
" * The LGPL. See the file COPYING.\n"
" */\n\n"
"#include \"defines.h\"\n\n";


static void
dump_array(float b[Nt][K])
{
  
  printf("static const float c2wideband_map[%d][%d] = {\n", Nt, K);
  int row, col;
  for (row = 0; row < Nt; row++ ) {
      printf("{ ");
      for (col = 0; col < K; col++ ) {
            printf("  %g", b[row][col]);
            if ( col < K - 1 )
                printf(", ");
            else
                printf(" }");
      }
      if ( row < Nt - 1 )
          printf(",\n");
      else
          printf("\n");
  }
  printf("};\n");
}



float
get_float(FILE * in, const char * name, char * * cursor, char * buffer,
 int size)
{
  for ( ; ; ) {
    char *	s = *cursor;
    char	c;

    while ( (c = *s) != '\0' && !isdigit(c) && c != '-' && c != '.' )
      s++;

    /* Comments start with "#" and continue to the end of the line. */
    if ( c != '\0' && c != '#' ) {
      char *	end = 0;
      float	f = 0;

      f = strtod(s, &end);

      if ( end != s )
        *cursor = end;
      return f;
    }

    if ( fgets(buffer, size, in) == NULL ) {
      fprintf(stderr, "%s: Format error. %s\n", name, format);
      exit(1);
    }
    *cursor = buffer;
  }
}

static void
load(FILE * file, const char * name, float b[Nt][K])
{
  char			line[1024];
  char *		cursor = line;

  *cursor = '\0';
  int row, col;

  for (row = 0; row < Nt; row++ ) {
      
      for (col = 0; col < K; col++ ) {
        
            b[row][col]  = get_float(file, name, &cursor, line, sizeof(line));
      }
  }
  
}


int main(int argc, char * * argv)
{
  float b[Nt][K];
  int	i;

  if ( argc < 2 ) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, format);
    exit(1);
  }

  i=0;
  FILE *	in = fopen(argv[i + 2], "r");

  if ( in == NULL ) {
      perror(argv[i + 2]);
      exit(1);
  }

  load(in, argv[i + 2], b);

  fclose(in);
 
  printf(header);
  
  printf("  /* %s */\n", argv[i + 2]);
  dump_array(b);

  return 0;
}
