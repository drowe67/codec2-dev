/* 
  FILE...: tollr.c
  AUTHOR.: David Rowe
  CREATED: July 2020

  Converts oneBitPerByte hard decisions to LLRs for LDPC testing.
*/

#include <stdio.h>
#include <stdint.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif /* _WIN32 */

int main(void) {
    uint8_t bit;
    
#ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
#endif /* _WIN32 */
    
    while(fread(&bit,sizeof(uint8_t), 1, stdin)) {
        float llr = 10.0*(1-2*bit);
        fwrite(&llr,sizeof(float),1,stdout);
    }
    return 0;
}
