function TransformedDomain(img::AbstractArray{<:ColorTypes.RGB, 2}, sigma_s::Float64, sigma_r::Float64)
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

    # Calculating the l1-normal distance of the pixels for each color channel
    # Σ|I'k(x)|
    # Where k is the color channel
    dx_normal_red[:,2:end] = dx_normal_red[:,2:end] .+ abs(dx_partial_red)
    dx_normal_green[:,2:end] = dx_normal_green[:,2:end] .+ abs(dx_partial_green)
    dx_normal_blue[:,2:end] = dx_normal_blue[:,2:end] .+ abs(dx_partial_blue)
    # Σ|I'k(Y)|
    dy_normal_red[2:end,:] = dy_normal_red[2:end,:] .+ abs(dy_partial_red)
    dy_normal_green[2:end,:] = dy_normal_green[2:end,:] .+ abs(dy_partial_green)
    dy_normal_blue[2:end,:] = dy_normal_blue[2:end,:] .+ abs(dy_partial_blue)

    # Calculating the horizontal and vertical derivatives of the transformed domain of each color channel
    # (1 + Σ|I'k(x)|) dx
    dx_horizontal_transf_red = zeros(size(red_channel))
    dy_vertical_transf_red = zeros(size(red_channel))
    dx_horizontal_transf_green = zeros(size(green_channel))
    dy_vertical_transf_green = zeros(size(green_channel))
    dx_horizontal_transf_blue = zeros(size(blue_channel))
    dy_vertical_transf_blue = zeros(size(blue_channel))
    # (1 + Σ|I'k(y)|) dy
    dx_horizontal_transf_red = ones(size(red_channel)) .+  (sigma_s ÷ sigma_r) * dx_normal_red
    dy_vertical_transf_red = ones(size(red_channel)) .+  (sigma_s ÷ sigma_r) * dy_normal_red
    dx_horizontal_transf_green = ones(size(green_channel)) .+  (sigma_s ÷ sigma_r) * dx_normal_green
    dy_vertical_transf_green = ones(size(green_channel)) .+  (sigma_s ÷ sigma_r) * dy_normal_green
    dx_horizontal_transf_blue = ones(size(blue_channel)) .+  (sigma_s ÷ sigma_r) * dx_normal_blue
    dy_vertical_transf_blue = ones(size(blue_channel)) .+  (sigma_s ÷ sigma_r) * dy_normal_blue

    # Returning the two derivatives
    return dx_horizontal_transf, dy_vertical_transf

end

function RecursiveFilter(img::AbstractArray{<:ColorTypes.RGB, 2}, red_derivative::Array{Float64, 2}, green_derivative::Array{Float64, 2}, blue_derivative::Array{Float64, 2},sigma_h::Float64)
    a = exp(-sqrt(2) ÷ sigma_h)
    ad_red = a.^red_derivative
    ad_green = a.^green_derivative
    ad_blue = a.^blue_derivative
    j::AbstractArray{<:ColorTypes.RGB, 2} = img

    # Getting the width of the image
    w = size(j, 2)

    # Equation 21 of the paper
    # J[n] = (1 - a^d)*I[n] + a^d*J[n - 1]
    for i in 2:w
        j.red[:, i] .= (1 .- ad_red[:, i]) .* j.red[:, i] + ad_red[:, i] .* (j.red[:, i-1] - j.red[:, i])
        j.green[:, i] .= (1 .- ad_green[:, i]) .* j.green[:, i] + ad_green[:, i] .* (j.green[:, i-1] - j.green[:, i])
        j.blue[:, i] .= (1 .- ad_blue[:, i]) .* j.blue[:, i] + ad_blue[:, i] .* (j.blue[:, i-1] - j.blue[:, i])
    end

    for i in w-1:-1:1
        j.red[:, i] .= (1 .- ad_red[:, i]) .* j.red[:, i] + ad_red[:, i] .* (j.red[:, i+1] - j.red[:, i])
        j.green[:, i] .= (1 .- ad_green[:, i]) .* j.green[:, i] + ad_green[:, i] .* (j.green[:, i+1] - j.green[:, i])
        j.blue[:, i] .= (1 .- ad_blue[:, i]) .* j.blue[:, i] + ad_blue[:, i] .* (j.blue[:, i+1] - j.blue[:, i])
    end

end

function TransposeImage(img::AbstractArray{<:Colors.RGB, 2})
    # Getting the height of the image
    h = size(img, 1)
    
    # Creating a transposed image with the correct size
    transposed_image = similar(img, size(img, 2), size(img, 1))
    
    for i in 1:h
        transposed_image[:, i] = img[i, :]
    end
    
    return transposed_image
end

using TestImages
using Images

img = testimage("mandrill")
display(img)

transposed_img = TransposeImage(img)
display(transposed_img)