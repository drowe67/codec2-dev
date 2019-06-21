/* debug_alloc.h
 *
 * Some macros which can report on malloc results.
 *
 * Enable with "-D DEBUG_MALLOC"
 */

#ifndef DEBUG_ALLOC_H
#define DEBUG_ALLOC_H

#include <stdio.h>

#ifdef DEBUG_ALLOC
// Debug calls

#ifdef CORTEX_M4
extern char * __heap_end;
register char * sp asm ("sp");
#endif

  static inline void * DEBUG_MALLOC(const char *func, size_t size) {
    void *ptr = malloc(size);
    fprintf(stderr, "MALLOC: %s %p %zd", func, ptr, size);
#ifdef CORTEX_M4

    fprintf(stderr, " : sp %p  heap_end %p", sp, __heap_end);
#endif
    if (!ptr) fprintf(stderr, " ** FAILED **");
    fprintf(stderr, "\n");
    return(ptr);
    }
  #define MALLOC(size) DEBUG_MALLOC(__func__, size)

  static inline void * DEBUG_CALLOC(const char *func, size_t nmemb, size_t size) {
    void *ptr = calloc(nmemb, size);
    fprintf(stderr, "CALLOC: %s %p %zd %zd", func, ptr, nmemb, size);
#ifdef CORTEX_M4
    fprintf(stderr, " : sp %p  heap_end %p", sp, __heap_end);
#endif
    if (!ptr) fprintf(stderr, " ** FAILED **");
    fprintf(stderr, "\n");
    return(ptr);
    }
  #define CALLOC(nmemb, size) DEBUG_CALLOC(__func__, nmemb, size)

  static inline void DEBUG_FREE(const char *func, void *ptr) {
    free(ptr);
    fprintf(stderr, "FREE: %s %p\n", func, ptr);
    }
  #define FREE(ptr) DEBUG_FREE(__func__, ptr)


#else //DEBUG_ALLOC
// Default to normal calls
  #define MALLOC(size) malloc(size)

  #define CALLOC(nmemb, size) calloc(nmemb, size)

  #define FREE(ptr) free(ptr)

#endif //DEBUG_ALLOC

#endif //DEBUG_ALLOC_H
