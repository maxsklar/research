A maintained version of this code (that will be further developed) is kept at:
https://github.com/maxsklar/BayesPy/tree/master/ConjugatePriorTools

This is simply a version that won't change for the purpose of the paper, in case the BayesPy version is changed and is no longer relevant to the paper.

Here's a command to try:
./sampleFromDirichletMultinomial.py -N 10000 -M 20 -A 1,5,2 | python findDirichletPrior.py -K 3 -M 100 -L INFO
cat test.csv | python findDirichletPrior.py