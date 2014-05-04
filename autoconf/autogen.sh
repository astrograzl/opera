#! /bin/bash
#
# A simple autogen
#
aclocal \
&& libtoolize \
&& automake --add-missing \
&& autoconf
