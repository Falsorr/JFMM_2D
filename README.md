# JFMM_2D
A implementation of the fast multipole method in 2D in Julia


# Usage of n.jl
The program uses a custom data Structure Point that has a "pos" field taking a complex number, it represents the position of the point in space, a "q" field representing its charge and a "pot" field representing the potential, this one should be initialised as 0 as the algorithmn aims to compute that value.
The program provides a function runFMM that takes as input a Vector of Point (data) and a Number for the desired precision (epsilon) that will modifiy each "pot" field of the Points in data.

# References
- https://math.nyu.edu/~greengar/shortcourse_fmm.pdf
- https://web.stanford.edu/class/cme324/classics/greengard-rokhlin.pdf
