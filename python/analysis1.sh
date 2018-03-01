#!/bin/sh
#BSUB -q 8nh
echo "Starting script."
NJOBS=100
source /afs/cern.ch/work/p/pvenkova/shipsetup1.sh

pathOut=/afs/cern.ch/work/p/pvenkova/disNoMagnetSBT/disNomagnetFiles
InputFile=ship.conical.muonDIS-TGeant4_${LSB_JOBINDEX}_rec.root 
#OutputFile=muonDISstudy_${LSB_JOBINDEX}.root
OutputFile=InCheck_${LSB_JOBINDEX}.root
GeoFile=geofile_full.conical.muonDIS-TGeant4.root 
File=/afs/cern.ch/work/p/pvenkova/disNoMagnetSBT/disNomagnetFiles/analysis_new.py
#analysNew.py
echo $pathOut
set -ux
echo "I am here"
if [ -f $pathOut/$OutputFile]; then
        echo "Target exists, nothing to do."
        exit 0
else
	echo "everything is correct , Now I am  executing  it"


	#xrdcp $pathOut/$File .
	xrdcp $pathOut/$InputFile .
	xrdcp $pathOut/$GeoFile .

	python $File   -f $InputFile -g $GeoFile 

	ls
	mv muonDISstudy.root muonDISstudy_${LSB_JOBINDEX}.root
	xrdcp $OutputFile $pathOut

	rm  $InputFile $OutputFile

fi 
