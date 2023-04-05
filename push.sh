#!/bin/bash

message="$1, `date`"
echo $message
git add --all .
git commit -m "$message"
git push origin master
