/*---------------------------------------------------------------------------*\

  FILE........: c2validate.h
  AUTHOR......: David Rowe
  DATE CREATED: 10 April 2013

  Encodes and decodes an array of speech samples using Codec 2 and compares
  it to a previously stored output to validate Codec operation.

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

#ifndef __C2VALIDATE__

int c2validate(int mode, short input_samples[], short output_samples[], char outfile[], int nsamples);

#endif
