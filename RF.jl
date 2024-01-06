a = [2 4; 6 16]

# Calculating the partial derivatives
dxP = diff(a, dims=2)
dyP = diff(a, dims=1)

# Creating the accumulated derivative matrixes
dxI = zeros(size(a))
dyI = zeros(size(a))

# Getting the size of the image matrix
rows = size(a, 1)
cols = size(a, 2)

# Calculating the l1-norm distance of the pixels
for i in 1:1
    for j in 2:cols
        for k in 1:rows
            dxI[k, j] = dxI[k, j] + abs(dxP[k, j-1])
        end
    end

    for j in 2:rows
        for k in 1:cols
            dyI[j, k] = dyI[j, k] + abs(dyP[j-1, k])
        end
    end
end

println("Original:")
for row in eachrow(a)
    println(row)
end

println("dxI:")
for row in eachrow(dxI)
    println(row)
end

println("dyI:")
for row in eachrow(dyI)
    println(row)
end