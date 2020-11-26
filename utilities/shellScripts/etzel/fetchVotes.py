#!/usr/bin/env python

import sys,urllib,re

instream = urllib.urlopen("http://www.doodle.com/wguhq4aq49utp439")
lines = instream.read()
instream.close()

dat=re.compile('td title=\"[^\"]+\"')
info=dat.findall(lines)

for x in info:
    print x
