#ifndef __CODEC2__
#define  __CODEC2__

#include <codec2/version.h>

#define CODEC2_MODE_3200 	0
#define CODEC2_MODE_2400 	1
#define CODEC2_MODE_1600 	2
#define CODEC2_MODE_1400 	3
#define CODEC2_MODE_1300 	4
#define CODEC2_MODE_1200 	5
#define CODEC2_MODE_700  	6
#define CODEC2_MODE_700B 	7
#define CODEC2_MODE_700C 	8
#define CODEC2_MODE_WB 		9
#define CODEC2_MODE_450 	10
#define CODEC2_MODE_450PWB 	11

struct CODEC2;

struct CODEC2 *  codec2_create(int mode);

#endif

