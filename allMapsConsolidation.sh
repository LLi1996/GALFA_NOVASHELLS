#!/bin/bash

# this shell script takes the _all_maps.pdf files from the individual 
# target folders and copy them to the /ALL folder where they'll be made into
# one single pdf


# information regarding the files
cycle="cycle1/"
fileEnding="_all_maps.pdf"
saveSuperDir="_all_targets/"
saveFolder="_fake_all_maps/"
saveName="fake_all_maps.pdf"

echo -e "\n\nThis script consolidates the files ending in" $fileEnding
echo -e "current cycle optoin:" $cycle

# saves the current path
superPath=$(pwd)
# all files are in the /Results dir
cd Results/

echo "Saving all files ending in" $fileEnding "to:"
echo -e "\t" $(pwd)/$saveSuperDir$saveFolder

cp fakeTargets/*/$cycle*$fileEnding $saveSuperDir$saveFolder

echo -e "\nNow we will merge all files ending in" $fileEnding "to one file"

cd $saveSuperDir
pdftk $saveFolder*$fileEnding cat output $saveName

echo -e "DONE\n"