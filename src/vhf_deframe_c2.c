
/*---------------------------------------------------------------------------*\

  FILE........: vhf_deframe_c2.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 8 March 2016

  C tool to extract codec2 data from freedv VHF 2400A/B/whatever frames
   

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2016 David Rowe

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
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "freedv_vhf_framing.h"

int main(int argc,char *argv[]){
    struct freedv_vhf_deframer * deframer;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    uint8_t *c2buf;
    uint8_t zbuf[] = {0,0,0,0,0,0,0};
    
    if(argc<3){
        fprintf(stderr,"usage: %s InputOneBitPerCharFile OutputC2File\n",argv[0]);
        exit(1);
    }
    
    /* Open files */
    if(strcmp(argv[1],"-")==0){
        fin = stdin;
    }else{
        fin = fopen(argv[1],"r");
    }
	
    if(strcmp(argv[2],"-")==0){
        fout = stdout;
    }else{
        fout = fopen(argv[2],"w");
    }

    /* Set up deframer */
    deframer = fvhff_create_deframer(FREEDV_VHF_FRAME_A,0);
    
    if(fin==NULL || fout==NULL || deframer==NULL){
        fprintf(stderr,"Couldn't open test vector files\n");
        goto cleanup;
    }
    
    /* allocate buffers for processing */
    bitbuf = (uint8_t*)malloc(sizeof(uint8_t)*96);
    c2buf = (uint8_t*)malloc(sizeof(uint8_t)*7);
    
    /* Deframe! */
    while( fread(bitbuf,sizeof(uint8_t),96,fin) == 96 ){
        if(fvhff_deframe_bits(deframer,c2buf,NULL,NULL,bitbuf))
            fwrite(c2buf,sizeof(uint8_t),7,fout);
        else
            fwrite(zbuf,sizeof(uint8_t),7,fout);
        
        if(fin == stdin || fout == stdin){
            fflush(fin);
            fflush(fout);
        }
    }
    
    free(bitbuf);
    free(c2buf);
    
    cleanup:
    fclose(fin);
    fclose(fout);
    fvhff_destroy_deframer(deframer);
    exit(0);
}

