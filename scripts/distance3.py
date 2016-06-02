#!/usr/bin/env python
# -*- coding=utf-8 -*-
import sys
import Levenshtein
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: You have to enter a source_word and a target_word'
        sys.exit(-1)
    source, target = sys.argv[1], sys.argv[2]
    print Levenshtein.distance(source, target)
