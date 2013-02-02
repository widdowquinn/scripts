#!/usr/bin/env python
#
# bisection.py
#
# Python module implementing the bisection method for root finding, based on
# OU M373 Unit I.1 and Numerical Methods in C Chapter 9.1
#
# The results of these two methods sometimes differ.  Interestingly, the more
# OU-like metod at http://deadline.3x.ro/bisection_method.html doesn't include
# the OU algorithm step (d), like the (interesting) Numerical Methods 
# algorithm.  Both proceed until the interval size is smaller than the 
# required accuracy which, as noted in the OU course does not guarantee that the
# rounded value of the root is correct to that accuracy - only that the 
# interval in which we found the root is smaller than the required accuracy.
# In some cases (such as Example 2.1), to be certain of the required result to 
# quoted accuracy we need to use the algorithm step (d) from the course text,
# - this indicates which of the two rounded options is correct.
#
# Possibly a blog post in this.
#
# (c) L.Pritchard 2012

from math import exp, log, cos, tan
import sys

def rtbis(fn, a, b, dp, iterlim=40, verbose=False):
    """ Using bisection, find the root of function fn known to lie between
        a and b. The root is refined until its accuracy is found to dp
        decimal places.  We run for a maximum of 40 iterations.

        This is the method described in Numerical Methods in C, Chapter 9.1,
        essentially transcribed into Python.
    """
    fa, fb = fn(a), fn(b)
    if fa * fb >= 0:
        raise ValueError, "Root must be bracketed for bisection"
    # We orient the search so that we're going 'uphill', with rtb as the
    # candidate root, and dx as the bisection interval; rtb is always negative
    rtb, dx = (a, b-a) if fa < 0 else (b, a-b)
    if verbose:
        print "\nNumerical Methods rtbis"
        hstr = "i(r)\trtb(ar)\txmid(cr)\tdx\tfn(xmid)"
        print hstr
        print '=' * len(hstr) + '=' * 4 * len(hstr.split('\t'))
    for i in range(iterlim):
        dx = dx * 0.5         # Halve the interval
        xmid = rtb + dx       # Find fn(midpoint)
        fmid = fn(xmid)
        if verbose:
            print "%d\t%f\t%f\t%f\t%f" % (i, rtb, xmid, dx, fmid)
        if fmid <= 0:         # If fn(midpoint) is negative, replace rtb
            rtb = xmid
        # If we've found the root, or if we have the desired accuracy, return
        # the root to the appropriate number of places
        if abs(dx) < 10**(-dp) or fmid == 0:
            root = round(rtb, dp)
            if verbose:
                print "Calculated root: %f" % rtb
                print "Root to %d d.p.: %s" % (dp, root)
            return root

def bisection(fn, a, b, dp, iterlim=40, verbose=False):
    """ Using the bisection algorithm described in M373 Unit I.1 Section 2,
        find the root of function fn known to lie between a and b.  The root
        is refined until its accuracy is found to dp decimal places.  We run
        for a maximum of 40 iterations.

        This method is written in slightly more idiomatic Python than rtbis
    """
    if verbose:
        print "\nOU Algorithm"
        hstr = "r\tar\tbr\tcr\tfar\tfbr\tfcr\tbr-ar"
        print hstr
        print '=' * len(hstr) + '=' * 6 * len(hstr.split('\t'))
    for r in range(iterlim):
        far, fbr = fn(a), fn(b)
        assert far * fbr < 0, "Root must be bracketed for bisection"
        # Algorithm steps (b)/(c)
        c = 0.5 * (a + b)
        fcr = fn(c)
        if verbose:
            print "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (r, a, b, c, 
                                                      far, fbr, fcr, b-a)
        if fcr == 0: # We have the root
            root = round(c, dp)
            break
        elif far * fcr < 0:
            b = c
        else:
            a = c
        # Algorithm step (d)
        if b - a <= 10**(-dp):
            a, b = round(a, dp), round(b, dp)
            if a == b:
                root = a
            else:
                c = 0.5 * (a + b)
                fcr = fn(c)
                if verbose:
                    print "%d*\t%f\t%f\t%f\t\t\t%f" % (r+1, a, b, c, fcr)
                if fcr == 0:
                    root = b
                elif fn(a) * fcr < 0:
                    root = a
                else:
                    root = b
            break
    if verbose:
        print "Root to %d d.p.: %s" % (dp, root)
    return root
    
            

if __name__ == '__main__':
    # Test function with OU exercises
    print "\n\nExample 2.1"
    def example_2_1(x):
        return x**3 - 0.75*x**2 - 4.5*x + 4.75
    rtbis(example_2_1, 1.5, 2, 3, verbose=True)
    bisection(example_2_1, 1.5, 2, 3, verbose=True)

    print "\n\nExercise 2.1"
    def exercise_2_1(x):
        return x*exp(-x) - 0.25
    rtbis(exercise_2_1, 0.3, 0.4, 3, verbose=True)
    bisection(exercise_2_1, 0.3, 0.4, 3, verbose=True)

    print "\n\nExercise 2.2"
    def exercise_2_2(x):
        return x * cos(x) - log(x)
    rtbis(exercise_2_2, 1, 1.6, 2, verbose=True)
    bisection(exercise_2_2, 1, 1.6, 2, verbose=True)

    print "\n\nExercise 2.3"
    def exercise_2_3(x):
        return x * tan(x) - log(x)
    rtbis(exercise_2_3, 3, 4, 5, verbose=True)
    bisection(exercise_2_3, 3, 4, 5, verbose=True)

    print "\n\nExercise 2.13"
    def exercise_2_13(x):
        return x**3 + 4*x**2 -10
    rtbis(exercise_2_13, 1, 2, 3, verbose=True)
    bisection(exercise_2_13, 1, 2, 3, verbose=True)

    print "\n\nExercise 2.14"
    def exercise_2_14(x):
        return x**3 - 5*x + 4
    rtbis(exercise_2_14, 1.25, 2, 5, verbose=True)
    bisection(exercise_2_14, 1.25, 2, 5, verbose=True)

    print "\n\nExercise 2.15"
    def exercise_2_15(x):
        return x * exp(-x) - 0.1
    rtbis(exercise_2_15, 3.31, 3.67, 5, verbose=True)
    bisection(exercise_2_15, 3.31, 3.67, 5, verbose=True)

