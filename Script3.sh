#!/bin/bash
echo "Enter the word you are looking for:"
read searchkey
find ./ -type f > filelist.txt
while read filename
do
wcount=$(grep $searchkey $filename | wc -l)
if [ $wcount != "0" ]
then
echo "*********** $filename **********"
grep $searchkey $filename
fi
done < filelist.txt
rm filelist.txt
