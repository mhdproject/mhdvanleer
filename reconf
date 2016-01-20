#!/bin/sh
rm -f config.cache
echo "- aclocal."
aclocal -I m4
echo "- autoconf."
autoconf
echo "- autoheader."
autoheader
echo "- automake."
automake -a
exit
