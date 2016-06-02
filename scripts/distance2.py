#!/usr/bin/env python
# -*- coding=utf-8 -*-

def lev_dist(source, target):
    if source == target:
        return 0

#words = open(test_file.txt,'r').read().split();

    # Prepare matrix
    slen, tlen = len(source), len(target)
    dist = [[0 for i in range(tlen+1)] for x in range(slen+1)]
    for i in xrange(slen+1):
        dist[i][0] = i
    for j in xrange(tlen+1):
        dist[0][j] = j

    # Counting distance
    for i in xrange(slen):
        for j in xrange(tlen):
            cost = 0 if source[i] == target[j] else 1
            dist[i+1][j+1] = min(
                            dist[i][j+1] + 1,   # deletion
                            dist[i+1][j] + 1,   # insertion
                            dist[i][j] + cost   # substitution
                        )
    return dist[-1][-1]

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print 'Usage: You have to enter a source_word and a target_word'
        sys.exit(-1)
    source, target = sys.argv[1], sys.argv[2]
    print lev_dist(source, target)
