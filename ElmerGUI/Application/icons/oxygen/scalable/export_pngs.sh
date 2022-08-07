#!/bin/bash

sizes="128 64 48 32 22 16"
folders="actions apps devices filesystems mimetypes"
smallexport="yes"
date=`date '+%F-%H-%M'`
curdir=$(pwd)'/'
icon=$1
usingGui="true"
max_small="48"

if [ "$icon" == "" ]; then
  icon=$(kdialog --getopenfilename $curdir)
else
  icon=$curdir$icon
  usingGui="false"
fi

iconName=$(basename $icon)
iconDir=$(basename `dirname $icon`)
iconPngName=$( echo $iconName | cut -d . -f -1 )".png"

inkscape --without-gui --export-png=$iconPngName --export-dpi=72 --export-background-opacity=0 --export-width=512 --export-height=512 $icon > /dev/null

for size in $sizes; do
  prefix="../${size}x${size}"
  # ====== shall we use a small icon if available?
  if [ $size -le $max_small ]; then
     smallicon="$iconDir/small/${size}x${size}/$iconName"
     if [ -e $smallicon ]; then
        inkscape --without-gui --export-png="../"${size}x${size}"/"$iconDir"/"$iconPngName --export-dpi=72 --export-background-opacity=0 --export-width=$size --export-height=$size $smallicon > /dev/null
     
     else
        convert -filter Sinc -resize ${size}x${size} $iconPngName "../"${size}x${size}"/"$iconDir"/"$iconPngName
     fi
  else
     convert -filter Sinc -resize ${size}x${size} $iconPngName "../"${size}x${size}"/"$iconDir"/"$iconPngName
  fi
  echo "Converted the icon named "$( echo $iconName | cut -d . -f -1 )" to size: " $size
done

rm $iconPngName

for size in $sizes; do
    svn add "../"${size}x${size}"/"$iconDir"/"$iconPngName
done

if $usingGui; then
  kdialog --msgbox "Icon converted and added to SVN"
else
  echo "Icon converted"
fi
