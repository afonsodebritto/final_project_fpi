using TestImages
img = testimage("mandrill")
using Images
display(img)

function transformed_domain(img::AbstractArray{<:ColorTypes.RGB, 2}, sigma_s::Float64, sigma_r::Float64)
    # Extracting each channel
    red_channel = red.(img)
    green_channel = green.(img)
    blue_channel = blue.(img)

    # Calculating the partial derivatives for each channel
    # I'k(x)
    # I'k(y)
    dx_partial_red = diff(red_channel, dims=2)
    dy_partial_red = diff(red_channel, dims=1)
    dx_partial_green = diff(green_channel, dims=2)
    dy_partial_green = diff(green_channel, dims=1)
    dx_partial_blue = diff(blue_channel, dims=2)
    dy_partial_blue = diff(blue_channel, dims=1)

    # Creating the accumulated derivative matrixes for each channel
    dx_normal_red = zeros(size(red_channel))
    dy_normal_red = zeros(size(red_channel))
    dx_normal_green = zeros(size(green_channel))
    dy_normal_green = zeros(size(green_channel))
    dx_normal_blue = zeros(size(blue_channel))
    dy_normal_blue = zeros(size(blue_channel))


    # Getting the size of the image matrix
    rows = size(img, 1)
    cols = size(img, 2)

    # Calculating the l1-normal distance of the pixels for each color channel
    # Σ|I'k(x)|
    # Σ|I'k(Y)|
    # Where k is the color channel
    for j in 2:cols
        for k in 1:rows
            dx_normal_red[k, j] = dx_normal_red[k, j] + abs(dx_partial_red[k, j-1])
            dx_normal_green[k, j] = dx_normal_green[k, j] + abs(dx_partial_green[k, j-1])
            dx_normal_blue[k, j] = dx_normal_blue[k, j] + abs(dx_partial_blue[k, j-1])
        end
    end

    dx_normal_red[:,2:end] = dx_normal_red[:,2:end] .+ abs(dx_partial_red)
    dx_normal_green[:,2:end] = dx_normal_green[:,2:end] .+ abs(dx_partial_green)
    dx_normal_blue[:,2:end] = dx_normal_blue[:,2:end] .+ abs(dx_partial_blue)

    for j in 2:rows
        for k in 1:cols
            dy_normal_red[j, k] = dy_normal_red[j, k] + abs(dy_partial_red[j-1, k])
            dy_normal_green[j, k] = dy_normal_green[j, k] + abs(dy_partial_green[j-1, k])
            dy_normal_blue[j, k] = dy_normal_blue[j, k] + abs(dy_partial_blue[j-1, k])
        end
    end

    # Calculating the scaling factor so it does not have to be calculated in every iteration
    scaling_factor::Float64 = sigma_s / sigma_r

    # Calculating the horizontal and vertical derivatives of the transformed domain of each color channel
    # (1 + Σ|I'k(x)|) dx
    # (1 + Σ|I'k(y)|) dy
    dx_horizontal_transf_red = zeros(size(red_channel))
    dy_vertical_transf_red = zeros(size(red_channel))
    dx_horizontal_transf_green = zeros(size(green_channel))
    dy_vertical_transf_green = zeros(size(green_channel))
    dx_horizontal_transf_blue = zeros(size(blue_channel))
    dy_vertical_transf_blue = zeros(size(blue_channel))

    dx_horizontal_transf_red = ones(size(red_channel)) .+  scaling_factor * dx_normal_red
    dy_vertical_transf_red = ones(size(red_channel)) .+  scaling_factor * dy_normal_red
    dx_horizontal_transf_green = ones(size(green_channel)) .+  scaling_factor * dx_normal_green
    dy_vertical_transf_green = ones(size(green_channel)) .+  scaling_factor * dy_normal_green
    dx_horizontal_transf_blue = ones(size(blue_channel)) .+  scaling_factor * dx_normal_blue
    dy_vertical_transf_blue = ones(size(blue_channel)) .+  scaling_factor * dy_normal_blue

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