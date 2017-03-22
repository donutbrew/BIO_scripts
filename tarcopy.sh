#!/bin/bash

dir=$1
dest=$2

tar -cf - $dir |tar -xf - -C $dest 2>&1
