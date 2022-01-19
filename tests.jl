include("n.jl")

#TESTS FOR N

# Ntests = [100, 1000, 5000, 10000, 40000, 80000, 100000]
Random.seed!(1)
# epsilonTest = 0.1
# pTest = Integer(ceil(log(sqrt(3), 1/epsilonTest)))


#io = open("testN.csv", "w");
#    write(io, "TESTS WITH N(epsilon = 0.1)\n")
#    write(io, "N, FMM_time (ns), BruteForceTime (ns)\n")
#
#    for i = 1:length(Ntests)
#        N = Ntests[i]
#        dataTest = []
#        for j = 1:N
#            x = rand(Float64)
#            y = rand(Float64)
#            z = rand(Float64)
#            push!(dataTest, PointFromCarthesian(x,y,z,1))
#        end
#        nTest = Integer(ceil(log(8, length(dataTest) / 4)))
#        if N < 10000
#            write(io, "$(length(dataTest)), $(runFMM(dataTest, nTest, pTest)), $(bruteForce(copy(dataTest)))\n")
#        else
#            write(io, "$(length(dataTest)), $(runFMM(dataTest, nTest, pTest))\n")
#        end
#    end
#close(io);

# TESTS FOR EPSILON

io = open("testE.csv", "w");
    NTest = 200
    EpsilonTests = Array(0.0001:0.001:0.1)
   write(io, "TESTS WITH EPSILON(N = $(NTest))\n")
   write(io, "Epsilon, FMM_time (ns)\n")
   d = []
   for j = 1:NTest
       x = rand(Float64)
       y = rand(Float64)
       z = rand(Float64)
       push!(d, PointFromCarthesian(x,y,z,1))
   end
   nTest = Integer(ceil(log(8, length(d) / 4)))

   for i = 1:length(EpsilonTests)
       e = EpsilonTests[i]
       pTest = Integer(ceil(log(sqrt(3), 1/e)))
        write(io, "$(e), $(runFMM(d, nTest, pTest))\n")
   end
close(io);

# d = []
# for j = 1:30000
#     x = rand(Float64)
#     y = rand(Float64)
#     z = rand(Float64)
#     push!(d, PointFromCarthesian(x,y,z,1))
# end
# println(bruteForce(d))
