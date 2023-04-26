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

    if [ ! -z "$2" ]
    then
        echo "Specify only the name of the build!"
    else
        mkdir $NAME
        mkdir $NAME/outputs

        cp ParameterFile.dat $NAME
        cp -r src/ $NAME/

        cp Makefile $NAME

        cd $NAME

        make

        echo "All done!"
    fi
fi