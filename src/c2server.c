/*---------------------------------------------------------------------------*\

  FILE........: c2server.c
  AUTHOR......: EA3IHI, David
  DATE CREATED: 28/05/2020

  

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#include "codec2.h"
#include "sine.h"
#include "dump.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <getopt.h>

#include <sys/select.h>
#include <sys/socket.h>

#include <netinet/in.h>
#include <netdb.h>


struct CODEC2 *codec2;
int serverPort = 1234;
int verbosity=0;

//Prints usage info.
void usage(char *argv0){
  printf("Usage: %s [OPTION]\n"
	 "\t-S          udp port to listen\n"
	 "\t-v          Verbose mode.\n",
	 argv0);
}


void c2Server(int portNumber) {
  struct sockaddr_in sa_read;
  //Initialize the codec2
  codec2 = codec2_create(CODEC2_MODE_2400);
  
  // open the socket
  int sockFd = socket(AF_INET, SOCK_DGRAM, 0);
  memset(&sa_read, 0x00, sizeof(struct sockaddr_in));
  sa_read.sin_family = AF_INET;
  sa_read.sin_port = htons(portNumber);
  sa_read.sin_addr.s_addr = htonl(INADDR_ANY);
  int ret = bind(sockFd, (struct sockaddr *)&sa_read, sizeof(struct sockaddr_in));

	if (verbosity) 
	{
        printf("Starting server %d\n", ret);
	}
	
	while (1) {
    unsigned char buffer[1024];

    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(sockFd, &fds);

    if (select(sockFd + 1, &fds, NULL, NULL, NULL) > 0) {
      if (verbosity)
        printf("data ready\n");
      socklen_t sl = sizeof(struct sockaddr_in);
      ssize_t n = recvfrom(sockFd, buffer, 1024, 0, (struct sockaddr *)&sa_read, &sl);
      if (verbosity)
        printf("got %d bytes\n", n);
      if (n==320) {
        unsigned char c2data[7];
        memset((unsigned char *)c2data, 0, sizeof(c2data));
		codec2_encode(codec2, c2data, (short *) buffer);
		c2data[6]= 0xC2;
		sendto(sockFd, c2data, sizeof(c2data), MSG_DONTWAIT, (struct sockaddr *)&sa_read, (ssize_t)sizeof(struct sockaddr_in));
        
		if (verbosity)
          printf("c2 sent\n");
      }
      if (n == 7) {
        short pcm[160];
        memset((unsigned char *)pcm, 0, sizeof(pcm));
        codec2_decode(codec2, (short *) pcm, (unsigned char*) buffer);
        sendto(sockFd, pcm, sizeof(pcm), MSG_DONTWAIT, (struct sockaddr *)&sa_read, (ssize_t)sizeof(struct sockaddr_in));
        if (verbosity)
          printf("pcm sent\n");
      }
    }
  }
}

int main(int argc, char **argv){
  int opt;
  char verb; //Main action of this run.
  
  verb = 'h';
  

  while((opt=getopt(argc,argv,"vS:"))!=-1){	
	printf("Processing flag: %c\n", opt);  
    switch(opt){
		case 'v':
			verbosity++;
			break;
		case 'S'://Socket
		  verb=opt;
		  serverPort = atoi(optarg);
		  if (serverPort == 0) {
			printf("Error, must specify port\n");
			exit(1);
		  }
			break;
	default:
		printf("Unknown flag: %c\n", opt);
		
	}
  }
  
  switch(verb){
	  case 'h'://Usage
		usage(argv[0]);
		exit(1);
	  case 'S':
		c2Server(serverPort);
		break;
	  default:
		printf("Usage error 2.\n");
		exit(1);
  }
	
  return 0;
}




