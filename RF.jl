sigma_s::float;
sigma_r::float;

a = [1 2 3;
     4 5 6;
     7 8 9]

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

# Calculating the scaling factor so it does not have to be calculated in every iteration
scaling_factor::float = sigma_s / sigma_r

# Calculating the horizontal and vertical derivatives of the transformed domain
dHx = zeros(size(a))
dVy = zeros(size(a))

for i in 1:rows
    for j in 1:cols
        dHx[i][j] = 1 +  (scaling_factor * dxI[i][j])
        dVy[i][j] = 1 +  (scaling_factor * dyI[i][j])
    end
end