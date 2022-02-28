echo "===== VQ10 codebook generation ========================="
echo "minimum maximum count"
# this creates minmax.csv file containing minimum and maximum values for LSPs, and stores quantization steps which default to 255
./src/c2sim ../VQ10/body.raw --lpc 10 --minmax 1 --verbose
echo "lsp vectors sample build"
# than we create lspvt.csv file where the dump of all frames is stored in a format of vectors
./src/c2sim ../VQ10/body.raw --lpc 10 --lspvt 1 --verbose
sudo mkdir csv
sudo rm ./csv/*
sudo cp minmax.csv ./csv
echo "looking for the best seeds"
# with --opti the lspvt.csv is read (the given number of vectors in --cdb) and than the distribution of vector occurences in analysed learn file is calculated
# the results are stored in lspvtopt.csv file from the most to the least frequent
# this file can be used as codebook however it is now treated as a entry point to LBG optimization
./src/c2sim ../VQ10/learn.raw --lpc 10 --cdb 100000 --opti 1 --verbose
sudo mv lspvt.csv ./csv
echo "test on the initial codebook"
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "first LBG pass"
# the --lbg option calculates the new centroids for each vector, on the end the set of centroids is written to the codebook.csv
# this file can be renamed to the lspvtopt.csv and works as better optimised codebook.
# the measure of codebook distortion is also calculated
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo mv lspvtopt.csv ./csv
sudo cp codebook.csv ./csv/codebook1.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "second LBG pass"
# subsequent iterations make optimisation better and better, showing the smaller codebook distortion
# when the drop in distortion is small the iterations can be stopped and resulting codebook assumed as optimal
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo cp codebook.csv ./csv/codebook2.csv
sudo rm lspvtopt.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "third LBG pass"
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo cp codebook.csv ./csv/codebook3.csv
sudo rm lspvtopt.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "fourth LBG pass"
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo cp codebook.csv ./csv/codebook4.csv
sudo rm lspvtopt.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "fifth LBG pass"
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo cp codebook.csv ./csv/codebook5.csv
sudo rm lspvtopt.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "sixth LBG pass"
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo cp codebook.csv ./csv/codebook6.csv
sudo rm lspvtopt.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "seventh LBG pass"
./src/c2sim ../VQ10/learn.raw --lpc 10 --vq10 1 --cdb 32768 --lbg --verbose
sudo cp codebook.csv ./csv/codebook7.csv
sudo rm lspvtopt.csv
sudo mv codebook.csv lspvtopt.csv
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "--------------------------------------------------------"
echo "final comparison:"
# now we will compare 3200, 1300 and VQ10 quantizers by given objective measures
# c2sim will provide SNR and SD measures
# embsd will provide EMBSD measure
echo "1. 3200 mode SQ"
./src/c2sim ../VQ10/test.raw --lpc 10 --lspd
./src/c2enc 3200 ../VQ10/test.raw test.c2
./src/c2dec 3200 test.c2 ./c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "2. 1300 mode SQ"
./src/c2sim ../VQ10/test.raw --lpc 10 --lsp
./src/c2enc 1300 ../VQ10/test.raw test.c2
./src/c2dec 1300 test.c2 ./c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
echo "3. VQ10 quantizer"
./src/c2sim ../VQ10/test.raw --lpc 10 --vq10 1 --cdb 32768 -o c2test.raw
./embsd ../VQ10/test.raw ./c2test.raw 1
sudo rm c2test.raw
sudo rm test.c2
echo "========================================================"
