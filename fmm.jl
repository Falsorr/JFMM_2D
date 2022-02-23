using Dates
using Random

mutable struct Point
    pos::Complex
    q::Number
    pot::Number
end

#Note here : c[l] = a_(l-1) in the local expansion, there is a change of variable
mutable struct Square
    center::Complex
    length::Number
    level::Integer
    parent
    children :: Vector{Square}
    points :: Vector{Point}
    ilist :: Vector{Square}
    nearPoints :: Vector{Point}
    Q::Number
    a::Vector{Number}
    b::Vector{Number}
end
function ComplexLog(z::Complex)::Complex
    #Purpose of this function is to compute log(z) with the main branch of the complex logarithmn function (for some reason it does not work in the code)
    p = log(z)
    while !(-pi < p.im <= pi)
        if p.im <= -pi 
            p += 1im*2pi
        else
            p += 1im*(-2pi)
        end
    end
    return p
end

function BaseSquare(center::Complex, length::Number)::Square
    return Square(center, length, 0, nothing, Vector{Square}(), Vector{Point}(),Vector{Square}(), Vector{Point}(), 0, Vector{Number}(), Vector{Number}())
end

function Split(SQUARES::Dict, square::Square)::Nothing
    newLength = square.length / 2
    newLevel = square.level + 1

    if !(newLevel in keys(SQUARES))
        SQUARES[newLevel] = Vector{Square}()
    end

    x = square.center.re
    y = square.center.im
    for i = -1:2:1
        for j = -1:2:1
            newCenter = x + i * newLength / 2 + (y + j * newLength / 2)*im
            child = Square(newCenter, newLength, newLevel, square, Vector{Square}(), Vector{Point}(), Vector{Square}(), Vector{Point}(), 0, Vector{Number}(), Vector{Number}())
            push!(square.children, child)
            push!(SQUARES[newLevel], child)
        end
    end
end

function PointInSquare(point::Point, square::Square)::Bool
    xC = square.center.re
    yC = square.center.im
    xP = point.pos.re
    yP = point.pos.im
    l = square.length / 2
    xCheck, yCheck = abs(xC - xP) <= l, abs(yC - yP) <= l

    return xCheck && yCheck
end

function SortSquares(squares::Matrix{Square}, square::Square, xmin::Integer, xmax::Integer, ymin::Integer, ymax::Integer)::Nothing
    # Function that sort the Squares in a 2D array so that we can access squares based on their position in space instead of their position in the data structure
    # SquareSorted[i][j] = i-th square on the x axis and j-th square on the y axis at the leaf level
    # Acknowledgment : Matas Keras helped greatly with this portion of the code by devising the function
    if (xmin == xmax) 
        squares[xmin, ymin] = square
        return
    end
    SortSquares(squares, square.children[1], xmin, Integer((xmax+xmin-1) / 2), ymin, Integer((ymin+ymax-1) / 2))
    SortSquares(squares, square.children[3], Integer((xmax+xmin+1) / 2), xmax, ymin, Integer((ymin+ymax-1) / 2))
    SortSquares(squares, square.children[2], xmin, Integer((xmax+xmin-1) / 2), Integer((ymin+ymax+1)/2), ymax)
    SortSquares(squares, square.children[4], Integer((xmax+xmin+1) / 2), xmax, Integer((ymin+ymax+1)/2), ymax)
end

function AssignPoints(SQUARES::Dict, data::Vector{Point}, n::Integer)::Nothing
    # We sort SQUARES so that we don't need to access each pair of point and squares
    # This reduces greatly time complexity compared at what was done before
    SquaresSorted = Matrix{Square}(undef, 2^n, 2^n)
    fill!(SquaresSorted, SQUARES[0][1])
    SortSquares(SquaresSorted, SQUARES[0][1], 1, 2^(n), 1, 2^(n))

    # Once the square are sorted, we iterate through each point and determine in which square they should belong 
    for i = 1:length(data)
        point = data[i]
        xCheck = floor(point.pos.re * 2^(n))+1
        yCheck = floor(point.pos.im * 2^(n))+1
        if xCheck == 2^(n) + 1
            xCheck = 2^(n)
        end
        if yCheck == 2^(n) + 1
            yCheck = 2^(n)
        end
        push!(SquaresSorted[Integer(xCheck), Integer(yCheck)].points, point)
        SquaresSorted[Integer(xCheck), Integer(yCheck)].Q += point.q
    end
    #We then transmit the information to the parents of each squares
    for i = n:-1:1
        squares = SQUARES[i]
        for j = 1:length(squares)
            square = squares[j]
            append!(square.parent.points, square.points)
            square.parent.Q += square.Q
        end
    end
end


function IsSeparated(square1::Square, square2::Square)::Bool
    isSperated = true
    v1 = GetVertices(square1)
    v2 = GetVertices(square2)

    for i =1:length(v1)
        if v1[i] in v2
            isSperated = false
            break
        end
    end

    return isSperated
end

function GetVertices(square :: Square)::Array
    vertices = []
    center = square.center
    l = square.length / 2

     for i = -1:2:1
        for j = -1:2:1
            push!(vertices, center + i * l +(j * l)im)
         end
     end
    return vertices
end

function GetIlistAndNearPoints(pSquaresSorted::Matrix{Square}, square::Square, plev::Integer)::Nothing
    #Find the parent in the sorted matrix
    pIndex = findall(x->x==square.parent, pSquaresSorted)[1]

    #We generate the submatrix of the near neighbours and second near neighbours of parent
    minXrange, maxXrange, minYrange, maxYrange = 2, 2, 2, 2
    if pIndex[1] == 1
        minXrange = 0
    elseif pIndex[1] == 2
        minXrange = 1
    end
    if pIndex[1] == 2^plev
        maxXrange = 0
    elseif pIndex[1] == 2^plev-1
        maxXrange = 1
    end
    if pIndex[2] == 1
        minYrange = 0
    elseif pIndex[2] == 2
        minYrange = 1
    end
    if pIndex[2] == 2^plev
        maxYrange = 0
    elseif pIndex[2] == 2^plev-1
        maxYrange = 1
    end

    #Submatrix containing parent's neighbour and second nearest neighbour
    NearParent = pSquaresSorted[pIndex[1]-minXrange:pIndex[1]+maxXrange, pIndex[2]-minYrange:pIndex[2]+maxYrange]

    #Once we have parent's near and second nearest neighbour, we find wich children are well separated of square (go into ilist), the ones that aren't form the nearPoints of square
    for i = 1:length(NearParent)
        for j = 1:4
            child = NearParent[i].children[j]
            dist = child.center - square.center
            if (-2*square.length <= dist.re <= 2*square.length) && (-2*square.length <= dist.im <= 2*square.length)
                append!(square.nearPoints, child.points)
            else
                push!(square.ilist, child)
            end
        end
    end
end


# a[l] = a_(l-1)
function MultipoleExpansion(square::Square, p::Integer)::Nothing
    #Compute mutipole expansion at mesh level
    points = square.points
    push!(square.a, 0 + 0im)
    square.a[1] += square.Q
    for k = 1:p
        push!(square.a, 0 + 0im)
        for i = 1:length(points)
            square.a[k+1] += -points[i].q * (points[i].pos - square.center)^k / k
        end
    end
end

function TranslationMultipoleExpansion(square::Square, p::Integer)::Nothing
    # Compute the translation of the multipole expansion of square's children to square
    children = square.children
    square.a = zeros(p+1)
    zM = square.center

    #Apply M2M operator to each children of square

    for i = 1:4
        child = children[i]
        zC = child.center
        z0 = zC-zM
        square.a[1] += child.a[1]
        for l = 1:p
            square.a[l+1] += - (child.a[1] * z0^l) / l
            for k = 1:l
                square.a[l+1] += child.a[k+1]*z0^(l-k) * binomial(l-1, k-1)
            end
        end
    end
end

#Note here : b[l] = b_(l-1) in the local expansion, there is a change of variable
function LocalExpansion(source::Square, target::Square, p::Integer)::Nothing
    a = source.a
    zS = source.center
    zT = target.center
    z = zS - zT
    target.b[1] += a[1] * ComplexLog(-z)
    for k = 1:p
        target.b[1] += (a[k+1])/(z^k) * (-1)^k
    end
    for l = 1:p
        target.b[l+1] += - a[1]/(l*z^l)
        for k = 1:p
            target.b[l+1] += 1/(z^l) * a[k+1] / (z^k) * binomial(l+k-1, k-1) * (-1)^k
        end
    end
end

#Note here : b[l] = a_(l-1) in the local expansion, there is a change of variable
function TranslationLocalExpansion(source::Square, target::Square)::Nothing
    zS = source.center
    zT = target.center
    z0 = zS - zT
    temp = copy(source.b)
    p = length(temp)-1
    if length(target.b) != length(temp)
        target.b = zeros(length(temp))
    end
    for j = 0:p-1
        for k = p-j-1:p-1
            temp[k+1] += - temp[k+2] * z0
        end
    end
    for j = 1:length(temp)
        target.b[j] += temp[j]
    end
end

function Init(center::Complex, l::Number, data::Vector{Point}, n::Number)::Dict
    println("Creating quadtree structure")
    SQUARES = Dict(0 => Vector{Square}())
    push!(SQUARES[0], BaseSquare(center, l))
    for i = 0:n-1
        squares = SQUARES[i]
        for j = 1:length(squares)
            Split(SQUARES, squares[j])
        end
    end
    println("Assigning Points")
    AssignPoints(SQUARES, data, n)
    println("Get Interactions Lists and Near Points")
    for i = 1:n
        # Sort the squares at the level above
        plev = i - 1
        pSquaresSorted = Matrix{Square}(undef, 2^plev, 2^plev)
        fill!(pSquaresSorted, SQUARES[0][1])
        SortSquares(pSquaresSorted, SQUARES[0][1], 1, 2^(plev), 1, 2^(plev))
        for j = 1:length(SQUARES[i])
            GetIlistAndNearPoints(pSquaresSorted, SQUARES[i][j], plev)
        end
    end
    return SQUARES
end

function UpwardPass(SQUARES::Dict, n::Integer, p::Integer)::Nothing
    println("Step 1 : Forming Multipole Expansions at mesh level")
    for i = 1:4^n
        MultipoleExpansion(SQUARES[n][i], p)
    end
    println("Step 2: Merging multipole expansions at coarser level using M2M operator")
    for l = n-1:-1:0
        for i = 1:4^l
            TranslationMultipoleExpansion(SQUARES[l][i], p)
        end
    end
end

function DownwardPass(SQUARES::Dict, n::Integer, p::Integer)::Nothing
    println("Step 3 : Forming local expansions at all but mesh level using M2L and L2L operators")
    #PROBABLY ERRORS HERE
    for l = 0:n-1
        squares = SQUARES[l]
        for i = 1:length(squares)
            square = squares[i]
            if length(square.b) == 0
                square.b = zeros(p+1)
            end
            for j = 1:length(square.ilist)
                otherSquare = square.ilist[j]
                LocalExpansion(otherSquare, square, p)
            end
            for j = 1:4
                child = square.children[j]
                if length(child.b) == 0
                    child.b = zeros(p+1)
                end
                TranslationLocalExpansion(square, child)
            end
        end
    end
    println("Step 4 : Forming local expansions at mesh level using M2L operator")
    #PROBABLY ERRORS HERE
    squares = SQUARES[n]
    for i = 1:length(squares)
        square = squares[i]
        for j = 1:length(square.ilist)
            otherSquare = square.ilist[j]
            LocalExpansion(square, otherSquare, p)
        end
    end
end

function EvaluationPotential_LocalExpansion(point::Point, square::Square, p::Integer)::Nothing
    zS = square.center
    b = square.b
    zP = point.pos
    for l = 1:p+1
        point.pot += (b[l] * (zP-zS)^(l-1)).re
    end
end

function Evaluationpotential_DirectEvaluation(point::Point, square::Square)::Nothing
    sourcePoints = square.nearPoints
    for i = 1:length(sourcePoints)
        if sourcePoints[i].pos != point.pos
            point.pot += (sourcePoints[i].q * (ComplexLog(point.pos - sourcePoints[i].pos))).re
        end
    end
end

function EvaluationPotential(squares::Vector{Square}, p::Integer)::Nothing
    #Uses the two  evaluation potential functions at each points, squares is SQUARES[n] (leaf level)
    println("Step 5 : Evaluating potential")
    for i = 1:length(squares)
        square = squares[i]
        for j = 1:length(square.points)
            point = square.points[j]
            EvaluationPotential_LocalExpansion(point, square, p)
            Evaluationpotential_DirectEvaluation(point, square)
        end
    end
end


#A function for running all steps of the FMM
function runFMM(data::Vector{Point}, epsilon::Number)::Nothing
    N = Integer(floor(log(4, length(data))))
    p = Integer(floor(log2(1/epsilon)))
    SQUARES = Init(0.5 + 0.5im, 1, data, N)
    UpwardPass(SQUARES, N, p)
    DownwardPass(SQUARES, N, p)
    EvaluationPotential(SQUARES[N], p)
end

#A function to compare to the bruteforce method
function bruteForce(data::Vector{Point})::Nothing
    for i = 1:length(data)
        for j = 1:length(data)
            if i != j
                data[i].pot += (data[j].q*log(data[i].pos-data[j].pos)).re
            end
        end
    end
end


#OUTDATED FUNCTIONS (NOT USED ANYMORE BECAUSE FOUND OPTIMIZED VERSION)
# function GetInteractionList(square::Square, SQUARES::Dict)::Nothing
#     # Get square's parent
#     parent = square.parent
#     pL = parent.level
#     parentList = Vector{Square}()
#     nearestParentList = Vector{Square}()
#     secondNearestParentList = Vector{Square}()
#     #Find nearest and second nearest neighbours of parent
#         #nearest
#         for i = 1:length(SQUARES[pL])
#             if !(IsSeparated(parent, SQUARES[pL][i]))
#                 push!(nearestParentList, SQUARES[pL][i])
#             end
#         end
#         #second nearest
#         for j = 1:length(nearestParentList)
#             nearNeighbour = nearestParentList[j]
#             for i = 1:length(SQUARES[pL])
#                 if !(IsSeparated(nearNeighbour, SQUARES[pL][i])) && nearNeighbour != parent
#                     push!(secondNearestParentList, SQUARES[pL][i])
#                 end
#             end
#         end
#         append!(parentList, append!(secondNearestParentList, nearestParentList))
#     #Add all children of squares in parentList that arent nearest or second nearest neighbours of square
#     potentialIList = Vector{Square}()
#     for i = 1:length(parentList)
#         for j = 1:length(parentList[i].children)
#             if !(parentList[i].children[j] in potentialIList) && parentList[i].children[j] != square
#                 push!(potentialIList, parentList[i].children[j])
#             end
#         end
#     end
#         #find nearest neighbours
#     nNeighbour = Vector{Square}()
#     for i = 1:length(potentialIList)
#         if !(IsSeparated(square, potentialIList[i])) && !(potentialIList[i] in nNeighbour) && potentialIList[i] != square
#             push!(nNeighbour, potentialIList[i])
#         end
#     end
#         #find points that aren't second nearest neighbours
#     for i = 1:length(potentialIList)
#         potSquare = potentialIList[i]
#         isNotClose = true
#         for j = 1:length(nNeighbour)
#             neighbour = nNeighbour[j]
#             if !IsSeparated(potSquare, neighbour)
#                 isNotClose = false
#             end
#         end
#         if isNotClose && !(potSquare in square.ilist)
#             push!(square.ilist, potSquare)
#         end
#     end
# end

# function GetNearPoints(square::Square, SQUARES::Dict)::Nothing
#     ilist = square.ilist
#     # Get square's neighbours and second neighbours
#     neighbours = Vector{Square}()
#     #Get the near neighbour
#     nearNeighbours = Vector{Square}()
#     secondNearNeighbours = Vector{Square}()
#     for i = 1:length(SQUARES[square.level])
#         otherSquare = SQUARES[square.level][i]
#         if !(IsSeparated(square, otherSquare)) && !(otherSquare in nearNeighbours)
#             push!(nearNeighbours, otherSquare)
#         end
#     end
#     append!(neighbours, nearNeighbours)
#     #second near neighbours
#     for i = 1:length(SQUARES[square.level])
#         otherSquare = SQUARES[square.level][i]
#         isClose = false
#         for j = 1:length(nearNeighbours)
#             neighbour = nearNeighbours[j]
#             if !(IsSeparated(otherSquare, neighbour))
#                 isClose = true
#             end
#         end
#         if isClose && !(otherSquare in nearNeighbours)
#             push!(secondNearNeighbours, otherSquare)
#         end
#     end
#     append!(neighbours, secondNearNeighbours)
#     #The neighbours and second neighbours are the squares not in the interaction list
#     for i = 1:length(neighbours)
#         otherSquare = neighbours[i]
#         if !(otherSquare in ilist)
#             append!(square.nearPoints, otherSquare.points)
#         end
#     end
# end
# function AssignPoints(SQUARES::Dict, data::Vector{Point}, n::Integer)::Nothing
#     squares = SQUARES[n]

#     for i = 1:length(data)
#         point = data[i]
#         for j = 1:length(squares)
#             square = squares[j]
#             if PointInSquare(point, square)
#                 push!(square.points, point)
#                 square.Q += point.q
#             end
#         end
#     end

#     for i = n:-1:1
#         squares = SQUARES[i]
#         for j = 1:length(squares)
#             square = squares[j]
#             append!(square.parent.points, square.points)
#             square.parent.Q += square.Q
#         end
#     end
# end
# function GetIlistAndNearPoints(SQUARES::Dict, square::Square)::Nothing
#     p = square.parent
#     #Get p's neighbour and second neighbours's children
#     plev = p.level
#     plen = p.length
#     slen = square.length
#     pC = p.center
#     psquares = SQUARES[plev]
#     pNeighboursChildren = Vector{Square}()
#     for i = 1:length(psquares)
#         potNeighbour = psquares[i]
#         dist = pC - potNeighbour.center
#         if (-2 * plen <= dist.im <= 2 * plen) && (-2 * plen <= dist.re <= 2 * plen)
#             append!(pNeighboursChildren, potNeighbour.children)
#         end
#     end

#     # If a square in pNeighboursChildren is a near or second nearest neighbour => points are added to the nearPoints
#     # if a square in pNeighboursChildren is not a near or second nearest neighbour => Added to ilist
#     for i = 1:length(pNeighboursChildren)
#         otherSquare = pNeighboursChildren[i]
#         dist = square.center - otherSquare.center
#         if (- 2 * slen <= dist.im <= 2 * slen) && (-2 * slen <= dist.re <= 2 * slen)
#             append!(square.nearPoints, otherSquare.points)
#         else
#             push!(square.ilist, otherSquare)
#         end
#     end
# end





