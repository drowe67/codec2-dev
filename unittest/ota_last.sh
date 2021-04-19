#!/bin/bash
# cat the log from from the latest test
cat `ls -td 2021* | head -n 1`/log.txt
