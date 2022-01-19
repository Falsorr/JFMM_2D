# Author: Daniel Fassler

using Dates
using LegendrePolynomials
using Random

mutable struct Point
    carthesian::Tuple
    spherical::Tuple
    charge::Number
    pot::Number
end

mutable struct Cube
    center :: Tuple
    length :: Number
    level :: Integer
    parent
    children
    points :: Array
    Q :: Number
    moments :: Dict
    hasBeenUsed :: Bool
    interactionList :: Array
    localExp :: Dict
    psihat :: Number
    nearPoints :: Array
end

function Distance(p1::Point, p2::Point)::Float64
    c1 = p1.carthesian
    c2 = p2.carthesian

    dist = sqrt((c1[1] - c2[1])^2+(c1[2] - c2[2])^2+(c1[2] - c2[2])^2)

    return dist
end

function PointFromCarthesian(x::Number, y::Number, z::Number, charge::Number)::Point
    cart = (x,y,z)
    if z != 0
        spher = (sqrt(x^2+y^2+z^2), atan(y,x), acos(z/sqrt(x^2+y^2+z^2)))
    else
        spher = (sqrt(x^2+y^2+z^2), atan(y,x), pi/2)
    end
    return Point(cart, spher, charge, 0)
end

function CreateBaseCube()
    return Cube((1/2,1/2,1/2), 1, 1, nothing, nothing, [], 0, Dict(), false, Vector{Cube}(), Dict(), 0, [])
end

function Split(cube::Cube)::Array
    length = cube.length / 2
    level = cube.level + 1
    x = cube.center[1]
    y = cube.center[2]
    z = cube.center[3]
    children = []
    for i = -1:2:1
        for j in -1:2:1
            for k in -1:2:1
                child = Cube((x + i*length/2, y + j*length/2, z + k * length/2), length, level, cube, nothing, [], 0, Dict(), false, Vector{Cube}(), Dict(), 0, [])
                push!(children, child)
            end
        end
    end
    cube.children = children
end

function IsSeparated(cube1::Cube, cube2::Cube)::Bool
    isSperated = true
    v1 = GetVertices(cube1)
    v2 = GetVertices(cube2)

    for i =1:length(v1)
        if v1[i] in v2
            isSperated = false
            break
        end
    end

    return isSperated
end

function PointInCube(point::Point, cube::Cube)::Bool
    isInCube = false

    xC, yC, zC = cube.center[1], cube.center[2], cube.center[3]
    xP, yP, zP = point.carthesian[1], point.carthesian[2], point.carthesian[3]
    l = cube.length / 2
    xCheck, yCheck, zCheck = abs(xC - xP) <= l, abs(yC - yP) <= l, abs(zC - zP) <= l

    if xCheck && yCheck && zCheck
        isInCube = true
    end
    return isInCube
end

LegendreFunction(n::Integer,m::Integer,x::Number)::Number = (-1)^m*(1-x^2)^(m/2) * dnPl(x, m, n)

SphericalHarmonics(n::Integer, m::Integer, theta::Number, phi::Number)::Complex = sqrt(factorial(big(n-abs(m))) / factorial(big(n+abs(m)))) * LegendreFunction(n, abs(m), cos(theta)) * exp(im*m*phi)

function MomentOfExpansion(n::Integer,m::Integer,points)::Complex
    moment = 0
    for i = 1:length(points)
        moment += points[i].charge * points[i].spherical[1] * SphericalHarmonics(n, -m, points[i].spherical[2], points[i].spherical[3])
    end
    return moment
end

Anm(n::Integer, m::Integer)::Number = (-1)^n / sqrt(factorial(big(n-m)) * factorial(big(n + m)))

function JmmPrime(m::Integer, mPrime::Integer)::Number
    if m*mPrime > 0
        return (-1)^mPrime * (-1)^(min(abs(mPrime), abs(m)))
    else
        return (-1)^mPrime
    end
end

function MomentsAtMesh(cubes::Array, p::Integer)::Nothing
    for i = 1:length(cubes)
        for j = 0:p
            for m = -j:1:j
                cubes[i].moments["n=$(j),m=$(m)"] = MomentOfExpansion(j, m, cubes[i].points)
            end
        end
    end
end


function GetVertices(cube :: Cube)::Array
    vertices = []
    center = cube.center
    l = cube.length / 2

     for i = -1:2:1
        for j = -1:2:1
             for k = -1:2:1
                 push!(vertices, center .+ (i * l, j * l, k * l))
             end
         end
     end

    return vertices
end


function getInteractionList(cube::Cube, CUBES::Dict)
    parent = cube.parent
    l = cube.level
    parentNeighbours = []
    for i = 1:length(CUBES["$(l-1)"])
        c = CUBES["$(l-1)"][i]
        if !(IsSeparated(parent, c))
            push!(parentNeighbours, c)
        end
    end
    interactionListCube = Vector{Cube}()
    for i = 1:length(parentNeighbours)
        neighbour = parentNeighbours[i]
        for j = 1:length(neighbour.children)
            cubeNeighbour = neighbour.children[j]
            if !(IsSeparated(cubeNeighbour, cube))
                push!(interactionListCube, cubeNeighbour)
            end
        end
    end
    return interactionListCube
end

function LocalExpansion(cube::Cube, p::Integer)::Dict
    localExp = Dict()
    x, y, z = cube.center[1], cube.center[2], cube.center[3]
    center = PointFromCarthesian(x, y, z, 1)
    rho, alpha, beta = center.spherical[1], center.spherical[2], center.spherical[3]
    moments = cube.moments
    for j = 0:p
        for k = -j:j
            localExp["j=$(j),k=$(k)"] = 0
            for i =1:p
                for m = -i:i
                    localExp["j=$(j),k=$(k)"] = localExp["j=$(j),k=$(k)"] + (moments["n=$(n=i),m=$(m)"]*JmmPrime(k, m)*Anm(i,m)*Anm(j,k)*SphericalHarmonics(j+n, m-k, alpha, beta))/(Anm(j+n, m-k)*rho^(j+n+1))
                end
            end
        end
    end
    return localExp
end

function computePsiHat(cube::Cube, p::Integer)::Number
    x,y,z,c = cube.center[1], cube.center[2], cube.center[3], 0
    center = PointFromCarthesian(x, y, z, c)
    r, theta, phi = center.spherical[1], center.spherical[2], center.spherical[3]
    psiHat = 0
    for j = 0:p
        for k = -j:j
            psiHat += cube.localExp["j=$(j),k=$(k)"] * SphericalHarmonics(j, k, theta, phi) * r^(j+1)
        end
    end
    return psiHat
end


#Form the interaction lists
function Init(data::Array, n::Integer)::Dict
    #Creating the Hierchical structure
    println("Initialization for FMM\n")
    println("Generating Octree structure")
    baseCube = CreateBaseCube()
    CUBES = Dict()
    CUBES["1"] = [baseCube]
    for l = 1:n-1
        children = []
        for i = 1:length(CUBES["$(l)"])
            cube = CUBES["$(l)"][i]
            append!(children, Split(cube))
        end
        CUBES["$(l+1)"] = children
    end

    # Assigning points and computing total charge
    println("Assigning points to structure")
    for l = n:-1:1
        cubes = CUBES["$(l)"]
        if l == n
            for i = 1:length(data)
                for j = 1:length(cubes)
                    if PointInCube(data[i], cubes[j])
                        push!(cubes[j].points, data[i])
                        cubes[j].Q += data[i].charge
                        break
                    end
                end
            end
        else
            for i = 1:length(cubes)
                for j = 1:length(cubes[i].children)
                    append!(cubes[i].points, cubes[i].children[j].points)
                    cubes[i].Q += cubes[i].children[j].Q
                end
            end
        end
    end
    println("Forming Interaction lists")
    for i = 3:n
        for j = 1:length(CUBES["$(i)"])
            cube = CUBES["$(i)"][j]
            cube.interactionList = getInteractionList(cube, CUBES)
        end
    end

    println("Forming Near Points\n")
    for i = 1:length(CUBES["$(n)"])
        cube = CUBES["$(n)"][i]
        for j = 1:length(CUBES["$(n)"])
            otherCube = CUBES["$(n)"][j]
            if !(IsSeparated(cube, otherCube))
                append!(cube.nearPoints, otherCube.points)
            end
        end
    end
    return CUBES
end

function FMM(CUBES::Dict,n ::Integer, p::Integer)::Nothing
    println("FMM\n")
    println("Upward Pass (Step 1/2)")
    # Upward pass
    for i = n:-1:1
        #Step 1 : Form Multipole Expansions at the finest level
        if i == n
            MomentsAtMesh(CUBES["$(i)"], p) # Compute moments of expansions at each level (except for 1 and 2)
        else
            #Step 2 : Merge Expansions at the coarser levels until the top level is reached
            for j = 1:length(CUBES["$(i)"])
                cube = CUBES["$(i)"][j]
                for k = 1:8
                    child = cube.children[k]
                    for o = 0:p
                        for m = -o:o
                            if "n=$(o),m=$(m)" in keys(cube.moments)
                                cube.moments["n=$(o),m=$(m)"] += child.moments["n=$(o),m=$(m)"]
                            else
                                cube.moments["n=$(o),m=$(m)"] = child.moments["n=$(o),m=$(m)"]
                            end
                        end
                    end
                end
            end
        end
    end

    # Downward pass
    println("Downard pass")
    println("Step 3")
    #Step 3 Form local expansions and adding them to psihat (all but leaf level)
    for l = 1:n-1
        for i = 1:length(CUBES["$(l)"])
            cube = CUBES["$(l)"][i]
            for m = 1:length(cube.interactionList)
                interactionCube = cube.interactionList[m]
                interactionCubeLocalExp = LocalExpansion(interactionCube, p)
                for j =0:p
                    for k = -j:j
                        if "j=$(j),k=$(k)" in keys(cube.localExp)
                            cube.localExp["j=$(j),k=$(k)"] += interactionCubeLocalExp["j=$(j),k=$(k)"]
                        else
                            cube.localExp["j=$(j),k=$(k)"] = interactionCubeLocalExp["j=$(j),k=$(k)"]
                        end
                    end
                end
                cube.psihat = computePsiHat(cube, p)
            end
        end

        for i = 1:length(CUBES["$(l)"])
            cube = CUBES["$(l)"][i]
            for j = 1:length(cube.children)
                child = cube.children[j]
                child.psihat += cube.psihat
            end
        end
    end

    println("Step 4")
    #Step 4 : Interactions at mesh level
    for i = 1:length(CUBES["$(n)"])
        cube = CUBES["$(n)"][i]
        for m = 1:length(cube.interactionList)
            interactionCube = cube.interactionList[m]
            interactionCubeLocalExp = LocalExpansion(interactionCube, p)
            for j =0:p
                for k = -j:j
                    if "j=$(j),k=$(k)" in keys(cube.localExp)
                        cube.localExp["j=$(j),k=$(k)"] += interactionCubeLocalExp["j=$(j),k=$(k)"]
                    else
                        cube.localExp["j=$(j),k=$(k)"] = interactionCubeLocalExp["j=$(j),k=$(k)"]
                    end
                end
            end
            computePsiHat(cube, p)
        end
    end

    println("Step 5/6")
    #Step 5/6 : Evaluate potential for particles inside the box and the particles inside the near neighbours
    for i = 1:length(CUBES["$(n)"])
        cube = CUBES["$(n)"][i]
        for j = 1:length(cube.points)
            point = cube.points[j]
            for k = 1:length(cube.nearPoints)
                if point != cube.nearPoints[j]
                    point.pot += 1 / Distance(point, cube.nearPoints[j])
                end
            end
        end
    end

    println("Step 7")
    #Step 7 : Add the far field potential to the near field potential
    for i = 1:length(CUBES["$(n)"])
        cube = CUBES["$(n)"][i]
        for j = 1:length(cube.points)
            point = cube.points[j]
            point.pot += cube.psihat
        end
    end
    println()
end

function runFMM(data::Array, n::Integer, p::Integer)::Nanosecond
    start = Dates.Time(Dates.now())
    CUBES = Init(data, n)
    FMM(CUBES, n, p)
    stop = Dates.Time(Dates.now())
    return stop - start
end

function bruteForce(data::Array)::Nanosecond
    println("BruteForce Method")
    start = Dates.Time(Dates.now())
    for i = 1:length(data)
        target = data[i]
        for j = 1:length(data)
            if j != i
                point = data[j]
                target.pot += 1/Distance(target, point)
            end
        end
    end
    stop = Dates.Time(Dates.now())
    println()
    return stop - start
end

if abspath(PROGRAM_FILE) == @__FILE__
    print("Running 'n.jl' as main file\n\n\n")
    println("Enter Accuracy")
    epsilon = parse(Float64, readline()) # Accuracy
    p = Integer(ceil(log(sqrt(3), 1/epsilon)))
    println("Enter Number of points")
    N = parse(UInt128, readline()) #Number of data points
    n = Integer(ceil(log(8, N/4))) #number of levels


    #Generating data points
    println("Generating Data")
    data = []
    for i = 1:N
        x = rand(Float64)
        y = rand(Float64)
        z = rand(Float64)
        c = rand(Float64)
        push!(data, PointFromCarthesian(x,y,z,c))
    end
            #FMM
    FMM_time = runFMM(data, n, p)

        #Brute force method
    dataCopy = copy(data)
    bruteForceTime = bruteForce(dataCopy)
    acc_test = true
    for i = 1:length(data)
        if (abs(data[i].pot - dataCopy[i].pot) > epsilon)
            acc_test = false
            break
        end
        println(data[i].pot,"\t", dataCopy[i].pot)
    end

    println("All potentials within Accuracy: $(acc_test)")
    println("FMM: It took $(FMM_time) to compute")
    println("Brute Force : It took $(bruteForceTime) to compute")
end
