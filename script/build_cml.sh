#!/bin/bash
wget http://www.iterativesolutions.com/user/image/cml.1.10.zip
unzip cml.1.10.zip
ls ${GITHUB_WORKSPACE}/octave/cml.patch
patch -p0 < ${GITHUB_WORKSPACE}/octave/cml.patch
cd cml/source
octave --no-gui -qf --eval "make"
