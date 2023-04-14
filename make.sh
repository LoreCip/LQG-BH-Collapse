#!/bin/bash

if [ "$1" == '-h' ]
then
    echo "Usage: ./make.sh OPT: path/build/folder"
else

    if [ -z "$1" ]
    then
        NAME=build
    else
        NAME=$1
    fi

    mkdir $NAME

    cp Makefile $NAME
    cp -r src/ $NAME/

    cd $NAME

    make

fi