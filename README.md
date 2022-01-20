# JFMM_3D
A implementation of the fast multipole method in 3D using spherical harmonics in Julia
THIS IS A WORK IN PROGRESS, THE PROGRAM ISN'T FINISHED AND SHOULDN'T BE TREATED AS A CORRECT IMPLEMENTATION 

# Usage of n.jl
The program in itself is a toy program that randomly generates data inside a cube of 1x1x1 but can be used to work with any dataset.
to use it, simply use the runFMM function by giving it as parameters :
- data : an 1D array containing Point objects, all points must be standardized to fit in a 1x1x1 cube
- p : an Integer limiting the number of terms inside the infinite sum
- n : the number of levels the FMM runs on
The function returns the time it took to do the FMM and modifies the pot values of the Point objects of data

# References
- https://math.nyu.edu/~greengar/shortcourse_fmm.pdf
- https://web.stanford.edu/class/cme324/classics/greengard-rokhlin.pdf
