function transformed_domain(img::Matrix{Float64}, sigma_s::Float64, sigma_r::Float64)

    # Calculating the partial derivatives
    # I'k(x)
    # I'k(y)
    dx_partial = diff(img, dims=2)
    dy_partial = diff(img, dims=1)

    # Creating the accumulated derivative matrixes
    dx_normal = zeros(size(img))
    dy_normal = zeros(size(img))

    # Getting the size of the image matrix
    rows = size(img, 1)
    cols = size(img, 2)

    # Calculating the l1-normal distance of the pixels for each color channel
    # Σ|I'k(x)|
    # Σ|I'k(Y)|
    # Where k is the color channel
    for i in 1:1
        for j in 2:cols
            for k in 1:rows
                dx_normal[k, j] = dx_normal[k, j] + abs(dx_partial[k, j-1])
            end
        end

        for j in 2:rows
            for k in 1:cols
                dy_normal[j, k] = dy_normal[j, k] + abs(dy_partial[j-1, k])
            end
        end
    end

    # Calculating the scaling factor so it does not have to be calculated in every iteration
    scaling_factor::Float64 = sigma_s / sigma_r

    # Calculating the horizontal and vertical derivatives of the transformed domain
    # (1 + Σ|I'k(x)|) dx
    # (1 + Σ|I'k(y)|) dy
    dx_horizontal_transf = zeros(size(img))
    dy_vertical_transf = zeros(size(img))

    for i in 1:rows
        for j in 1:cols
            dx_horizontal_transf[i, j] = 1 +  (scaling_factor * dx_normal[i, j])
            dy_vertical_transf[i, j] = 1 +  (scaling_factor * dy_normal[i, j])
        end
    end

    # Returning the two derivatives
    return dx_horizontal_transf, dy_vertical_transf

end

a = rand(Float64, 3, 3)  # Example input matrix
println("Original Matrix:")
for row in eachrow(a)
    println(row)
end

result_dx, result_dy = transformed_domain(a, 1.0, 1.0)

println("\nResult for dx_horizontal_transf:")
for row in eachrow(result_dx)
    println(row)
end

println("\nResult for dy_vertical_transf:")
for row in eachrow(result_dy)
    println(row)
end