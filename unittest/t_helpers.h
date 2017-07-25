/* 
 * File:   t_helpers.h
 * Author: phil
 *
 * Created on 21 July 2017, 14:20
 */

#ifndef T_HELPERS_H
#define	T_HELPERS_H

void test(char * tfn);
void test_failed();
void test_failed_s(char * expected, char * res);
void test_failed_f(float expected, float res);

char *fn;


#endif	/* T_HELPERS_H */

