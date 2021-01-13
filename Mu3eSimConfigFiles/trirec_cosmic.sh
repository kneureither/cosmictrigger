#!/bin/sh

if [ -x "./mu3eTrirec/mu3eTrirec" ] ; then
    trirec="./mu3eTrirec/mu3eTrirec"
elif [ -x "./trirec" ] ; then
    trirec="./trirec"
fi

"$trirec" --cosmic --conf="trirec_cosmic.conf" "$@"
