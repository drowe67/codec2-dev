#! /usr/bin/python3

""" Find stack usage using
 - compiler generated tables of stack use per function, *.c.su
 - run time trace output (function_trace.out) from compiler added enter/exit calls.

Just for ofdm_stack at this point
"""

COMP_DIR = 'unittest/CMakeFiles/ofdm_stack.dir'
EXE_FILE = 'unittest/ofdm_stack'

import sys
import pathlib
import subprocess

##########################
# Data Structures
su_data = {} # <function> : <stack_size>

##########################
# Read compiler generated tables of stack use per function, *.c.su
#
#  ofdm_stack.c:184:6:dummy_code	16	static
#

p = pathlib.Path(COMP_DIR)
for fpath in p.glob('**/*.c.su'):
    with fpath.open() as f:
        for line in f.readlines():
            try:
                words = line.split()
                size = int(words[1])
                words = words[0].split(':')
                su_data[words[3]] = size
            except: pass # skip this line if there are errors

##########################
# Read trace file, convert addresses to names, track stack

max_stack_depth = 0
cur_stack_depth = 0
stack = []  # List of tuples of (function names, cur_stack_depth)
last_func = 'start'

def walk_stack():
    trace = ''
    for entry in stack:
        trace += entry[0] + ' '
    return(trace)

# Open trace
with open("function_trace.out", "r") as f:
    for line in f.readlines():
        #print('Line: "{}"'.format(line.strip()))
        words = line.split()
        # Note addr2line needs addr in hex!
        addr = words[1]
        # ignore system calls
        if (int(addr, 0) < 0x7f00000000):
            if (words[0] == 'e'):
                result = subprocess.run(['addr2line', '-f', addr, '-e', EXE_FILE],
                                    stdout=subprocess.PIPE)
                result.check_returncode()
                # function name is first line of stdout
                if (result.stdout):
                    lines = result.stdout.decode().split('\n')
                    func = lines[0].strip()
                else: sys.error('unknown function at address {}'.format(addr))

                # Push last info
                stack.append((last_func, cur_stack_depth))
                last_func = func

                # Update
                cur_stack_depth += su_data[func]
                #print('func: "{}" = {}'.format(func, cur_stack_depth))
                if (cur_stack_depth > max_stack_depth):
                    max_stack_depth = cur_stack_depth
                    max_stack_trace = walk_stack()

                # end if ('e')
            elif (words[0] == 'x'):
                # Pop
                (last_func, cur_stack_depth) = stack.pop()
                #print('pop:  "{}" = {}'.format(last_func, cur_stack_depth))

print('Max Stack Depth = {}'.format(max_stack_depth))
print('Max Stack at: {}'.format(max_stack_trace))
