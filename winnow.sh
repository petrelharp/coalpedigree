#/bin/bash

if [ -z $1 ]
then
    cat <<-"EOF"
        winnow.sh: restrict to blocks above a certain length.

        Usage:
            ./winnow.sh (minlength) (infile|-)
        where (minlength) is a number.

EOF
    exit 0
fi

read HEADER
echo $HEADER
awk '{ if ( $4-$3 > '$1' ) print }'
