# 2014-Math4620-Project2
Crout factorization vs. Gauss-Seidel iteration on a particular tridiagonal matrix (Homework for Numerical Analysis 2)

- Function pointer notation in C finally clicked for me
- I basically implemented a class with a little polymorphism in C, which seems silly now.
- I had fun optimizing this
  - cutting 5 hours down to under a minute was a bigger improvement than expected
  - unrolling the first and last few iterations to avoid edge-testing conditionals inside the loop was probably overkill
- There are a few other fun/silly things, like the data structures of sparse matrices and `crash()`
