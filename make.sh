#!/bin/bash

if [ "$1" == '-h' ]
then
    echo "Usage: ./make.sh [remove path/to/build/folder] OPT: path/to/custom/build/folder"

elif [ "$1" == 'remove' ]
then

    DIRECTORY=$2

    if [ ! -d "$DIRECTORY" ]
    then
        echo "'$DIRECTORY' does not exist. Cannot remove."
    else
        rm -rf $2
    fi
else

    if [ -z "$1" ]
    then
        NAME=build
    else
        NAME=$1
    fi

    if [ ! -z "$3" ]
    then
        echo "Specify only the name of the build!"
    else
        mkdir -p $NAME
        mkdir -p $NAME/outputs

        cp ParameterFile.par $NAME
        cp -r src/ $NAME/

        cp Makefile $NAME

        cd $NAME

        make  ## NOHDF5=1

        echo "All done!"
    fi
fi