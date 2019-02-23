#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>

#include "semihosting.h"

extern void initialise_monitor_handles(void);
extern int errno;

int semihosting_init(void) {

    /* ARM Standard semihosting, newlib implementation.
     *
     * Without C startup we need to call initialise_monitor_handles() here.
     *
     * This code opens stdout and stderr to ":tt" which has special meaning
     * to the Arm debugger, use the debugger's text console.
     *
     * The st-util implementation does not recognize that, so it just trys to
     * open two files pointing at ":tt" which does not work.
         *  
     * So reopen them here to 2 separate files.  (Which may be named pipes
     * depending how this is being run from the PC.)
     *
     * Freopen will call fclose -> fflush -> lseek if the file is buffered.
     * This will result in some semihosting calls that are not implemented
     * in st-util there will be error messages.
     * To avoid those, set the files to unbuffered.
     */

    initialise_monitor_handles();    

    stdout = freopen("stm_stdout.txt", "w", stdout);
        if (!stdout) {
            fprintf(stderr, "Error %d reopening stdout\n", errno);
            return(errno);
        }
    setbuf(stdout, NULL);

    stderr = freopen("stm_stderr.txt", "w", stderr);
        if (!stderr) {
            fprintf(stdout, "Error %d reopening stderr\n", errno);
            return(errno);
        }
    setbuf(stderr, NULL);

    return(0);

}

/* vi:set ts=4 et sts=4: */
