#!/bin/bash

## script to format the source code uniformly.
## Only works under linux. If you are using OSX, replace
## 'indent -gnu' with 'gnuindent'
indent -gnu -npsl -npcs -nbs -nsaf -nsai -nsaw -nprs -bap -pmt -nut -lp -hnl -nut -l200 ./src/*/*.c
indent -gnu -npsl -npcs -nbs -nsaf -nsai -nsaw -nprs -bap -pmt -nut -lp -hnl -nut -l200 ./src/*/*/*.c
