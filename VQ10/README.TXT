# codec2_VQ10
Experiment tool for new Vector Quantizer for use in Codec2, 

this version 2.0 implements several changes to the start version:
- new options and functionality of c2sim
- new learning speech data  prepared from TSP Lab samples
- new EMBSD quality measurement for additional check
- Linde-Buzo-Gray algorithm for optimising the codebook

NOTES:
- from now lspvt.csv, lspvtopt.csv and codebook.csv files share the same format
- you do not put any data to the files manually, it is provided by options
- the changes has been made to c2sim.c, defines.h, phase.c, sine.c, sine.h

to get new environment you need to rebuild the c2sim:

cd codec2/build_linux
make c2sim

than you will need to compile embsd

gcc -o embsd embsd.c -lm

you should have codec2/VQ10 directory and there the raw speech samples for this branch are stored. we will use them
to start new version of VQ10 you will perform a build of it using c2sim.
to be sure youre going from scratch please

sudo rm *.csv

please prepare csv directory 

mkdir csv

than you can run script which will perform all steps for building the vq10 codebook
I suggest to run it in background as it will take some time (for me it is about 1-2 hrs.)

sudo ./lbg.sh &> out.txt &

to have and idea what is going on you can inspect script,  it is damn simple and comments show the idea
finally you can listen to the test sample. it is better now because thaere are a lot of short frazes of many speakers in a short file

list the results

cat out.txt

and read the objective parameters
subsequent codebooks are stored in ./csv directory

this will give you idea of how it sounds originally:

aplay -f S16_LE ../VQ10/test.raw

than try 3200 mode quantization (but without decimation):

./src/c2sim ../VQ10/test.raw --lpc 10 --lspd -o - | aplay -f S16_LE

than try 1300 mode quantization (but without decimation):

./src/c2sim ../VQ10/test.raw --lpc 10 --lsp -o - | aplay -f S16_LE

and finally try VQ10 quantization (without decimation as well):

./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o - | aplay -f S16_LE

this version of VQ10 needs 15 bits to code the LSPs for the frame.
for --vq10 1 option to run properly there are two files needed:
minmax.csv (minimum and maximum values for LSP rescaling)
lspvtopt.csv (the codebook)

regarding c2sim extensions what might be interesting are the new parametrization to play with:
--trsh 10.0    sets the voiced/unvoiced treshold to the given number
--trn          will perform voicing phases for lower half of harmonics and unvoiced for upper half
--nw 220       sets the analyse window width (default is 270) 
