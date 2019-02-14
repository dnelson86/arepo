#!/bin/bash


DATE=$(git log -n 1 2> /dev/null | head -n4 | grep "Date" | cut -d' ' -f4-)
COMMIT=$(git log -n 1 2> /dev/null | head -n1 | grep "commit" | cut -d' ' -f2-)


if  [[ $DATE == "" ]]
    then
    DATE="unknown"
fi

if  [[ $COMMIT == "" ]]
    then
    COMMIT="unknown"
fi



if [ -f $BUILD_DIR/version.c ]
    then

    COMMIT2=$(grep "GIT_COMMIT" $BUILD_DIR/version.c)

    if [[ $COMMIT2 == *$COMMIT* ]] #it's the same commit
    then
        exit
    fi
fi

cp $SRC_DIR/gitversion/version $BUILD_DIR/version.c
sed -i.bu "s/_DATE_/$DATE/g" $BUILD_DIR/version.c
sed -i.bu "s/_COMMIT_/$COMMIT/g" $BUILD_DIR/version.c
rm $BUILD_DIR/version.c.bu
