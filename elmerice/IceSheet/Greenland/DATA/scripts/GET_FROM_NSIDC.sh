#!/bin/bash

## get data from nsidc
if [ $# -eq 0 ]
  then
    echo "No arguments supplied -- ABORT --"
    echo
    echo "provide the url of the file to download"
fi
## test if authentification file exist
if [ ! -e "$HOME/.netrc" ]; then
  echo "~/.netrc does not exist"
  echo 
  echo "see https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget"
  echo 
  echo "  echo "machine urs.earthdata.nasa.gov login \<uid\> password \<password\>" >> ~/.netrc"
  echo "  chmod 0600 ~/.netrc"
  echo 
  echo "or directly download file from you web browser from this link:"
  echo $1 
  echo 
  exit 1
fi
# create cookies file
if [ ! -f "~/.urs_cookies" ]; then
   touch ~/.urs_cookies
fi
## download file using curl
curl -O -b ~/.urs_cookies -c ~/.urs_cookies -L -n $1


