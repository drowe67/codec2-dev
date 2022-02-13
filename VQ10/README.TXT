# codec2_VQ10
Experiment tool for new Vector Quantizer for use in Codec2, this version focused on 10 dim VQ

all changes were done only to c2sim.c

so to get new environment just:

cd codec2/build_linux
make c2sim

you should have codec2/VQ10 directory and there the raw speech samples for this branch are stored. we will use them
to start new version of VQ10 you will perform a build of it using c2sim.
to be sure youre going from scratch please

sudo rm *.csv

first step will be to create minmax.csv storing first step quantizer resolution and minimum and maximum values for each LSP in the sample.
just run

./src/c2sim ../VQ10/body.raw --lpc 10 --minmax 1 --sqdd 255

this will create minmax.csv file you can inspect by

sudo nano ./minmax.csv

first line will hold resolution given in sqdd option (if you omit this it will be 255 by default.)
and than you will have 10 lines of two float values separated by colons.
ok so far so good. now you do the vector database by recording all vectors of given sample.
(I suggest to work on the same sample for minmax and lspvt however it is not an obligation).
to do this please run:

./src/c2sim ../VQ10/body.raw --lpc 10 --lspvt 1

you should get lspvt.csv file where all vectors of a sample will be written in CSV format.
the lspvt value has to be done. this is the version of current "event definition" and "error calculation".
please just put 1 the others are under dev
now the stupid part. sorry for that. quick and dirty, it was just quicker to do for me.
you have to read the number of frames the preavious command will show and put the lower number in the first line of the resulting file
I normally do it by

sudo nano ./lspvt.csv

and manually add one number. in the case of body.raw where there is about 99k lines I put 90000 or 95000.
for longer samples please dont put a number bigger than 100000.
now as you did all this is the most nice part. we will learn or rather sort vectors calculating the distribution of closest vectors in another sample.
to do this you do:

./src/c2sim ../VQ10/learn.raw --lpc 10 --opti 1

make a good coffee or a good tee, in my case Twinnings Assam will do the job. it will take time of this cup of tee.
the opti value has to be the same 1. but this command reads minmax.csv and lspvt.csv than for each vector in learn.raw it will find a closest vector in the database of vectors in lspvt.csv limited to the number in the first line. the command will create lspvtopt.csv file where you will have another table.
please note the number of "maenigfull vectors" as an output of the command. 
the table consists of a vectors but in the first position you have the number of hits in learn.raw sample.so you will see the distribution too.

and now you are finally ready to hear the vq10 quantizer. please take:

./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32000 -o - | aplay -f S16_LE

you should hear the voice transcoded by VQ10. on the end you will get SD measurement too. you have to put this silly 1 too.
the cdb option is the number of vectors read from the lspvtopt.csv into memory to search for the closest.
but you can run all other old options too so 

./scr/c2sim ../VQ10/test.raw --lpc 10 --lsp -o - | aplay -f S16_LE

will work too as ALL other old working options of c2sim.

I published some experiment results in the codec2 repository discussion.
If you want to store the VQ10 "version" you have to make some new directory and copy all created CSV files into it.
Than you can start with a new combination of parameters like first step resolution, body and learned vector databases.
For the larger files of body I just split the result in lspvt.csv into 100000 lines files using split.
Than I optimize them one by one. There is no practical limit of learn sample however bigger files will take some more time to calculate.
The trends of experiment published on the same discussion

