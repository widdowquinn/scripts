# Short, ad hoc script to find a solution to a puzzle 
# from Douglas Hofstadter's book, Metamagical Themas
#
# L. Pritchard 2013

import itertools
import random

def find_hofstadter(l):
    l0 = None
    while l != l0:
        l0 = l[:]
        digits = [int(j) for j in \
                      list(itertools.chain.from_iterable([list(str(e)) \
                                                              for e in l]))]
        for i in range(len(l)):
            l[i] = digits.count(i) + 1
    return l

for i in range(10):
    print i, find_hofstadter([i]*10)

#for i in range(20):
#    l = random.sample(range(10), 10)
#    print l, find_hofstadter(l)

print find_hofstadter([4, 10, 10, 2, 9, 5, 1, 7, 3, 1])
