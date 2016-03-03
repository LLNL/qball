#!/usr/bin/env python


from bs4 import BeautifulSoup
import sys

energy_tol = 1e-5;
verbose = False;

#xml parsing
filename1 = sys.argv[1];
filename2 = sys.argv[2];
file1 = BeautifulSoup(open(filename1, 'r').read(), "lxml");
file2 = BeautifulSoup(open(filename2, 'r').read(), "lxml");

iterations1 = file1.findAll("iteration");
iterations2 = file2.findAll("iteration");

correct = True;

if(len(iterations1) != len(iterations2)):
    print "Error: different number of iterations";
    exit(1);

for index, (it1, it2) in enumerate(zip(iterations1, iterations2)):

    iteration_correct = True;

    print "Iteration", index, "\t:\t",

    etotal1 = float(it1.etotal.text);
    etotal2 = float(it2.etotal.text);

    if(abs(etotal1 - etotal2) > energy_tol): iteration_correct = False;

    if(verbose): print "etotal", etotal1 - etotal2;

    if(iteration_correct):
        print "[   OK   ]";
    else:
        print "[  FAIL  ]";
    
    correct = correct and iteration_correct

if (not correct): exit(1);
