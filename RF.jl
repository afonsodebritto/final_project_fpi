function TransformedDomain(img::AbstractArray{<:ColorTypes.RGB, 2}, sigma_s::Float32, sigma_r::Float32)
    # Calculating the partial derivatives for all the channels
    # I'k(x)
    # I'k(y)
    I = float.(img)
    dx_partial = diff(I, dims=2)
    dy_partial = diff(I, dims=1)

    # Creating the accumulated derivative matrixes for all the channels
    h = size(I, 1)
    w = size(I, 2)
    dx_normal::Matrix{Float32} = zeros(Float32, h, w)
    dy_normal::Matrix{Float32} = zeros(Float32, h, w)

    # Calculating the l1-normal distance of the pixels
    # Σ|I'k(x)|
    # Where k is the color channel
    dx_partial = float.(channelview(dx_partial))
    for i in 1:3
        dx_normal[:, 2:end] = dx_normal[:, 2:end] + abs.(dx_partial[i, :, :])
    end

    # Σ|I'k(Y)|
    # Where k is the color channel
    dy_partial = float.(channelview(dy_partial))
    for i in 1:3
        dy_normal[2:end, :] = dy_normal[2:end, :] + abs.(dy_partial[i, :, :])
    end

    # Calculating the horizontal and vertical derivatives of the transformed domain
    # (1 + (σs÷σr)*Σ|I'k(x)|) dx
    dx_horizontal::Matrix{Float32} = ones(Float32, h, w)
    sigma_s = 0
    sigma_r = 1

    dx_horizontal += (sigma_s ÷ sigma_r) * dx_normal

    # (1 + Σ|I'k(y)|) dy
    dy_horizontal::Matrix{Float32} = ones(Float32, h, w)
    dy_horizontal += (sigma_s ÷ sigma_r) * dy_normal

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
        red.(j[:, i]) .= (1 .- ad_red[:, i]) .* red.(j[:, i]) + ad_red[:, i] .* (red.(j[:, i-1]) - red.(j[:, i]))
        green.(j[:, i]) .= (1 .- ad_green[:, i]) .* green.(j[:, i]) + ad_green[:, i] .* (green.(j[:, i-1]) - green.(j[:, i]))
        blue.(j[:, i]) .= (1 .- ad_blue[:, i]) .* blue.(j[:, i]) + ad_blue[:, i] .* (blue.(j[:, i-1]) - blue.j([:, i]))
    end

    for i in w-1:-1:1
        red.(j[:, i]) .= (1 .- ad_red[:, i]) .* red.(j[:, i]) + ad_red[:, i] .* (red.(j[:, i+1]) - red.(j[:, i]))
        green.(j[:, i]) .= (1 .- ad_green[:, i]) .* green.(j[:, i]) + ad_green[:, i] .* (green.(j[:, i+1]) - green.(j[:, i]))
        blue.(j[:, i]) .= (1 .- ad_blue[:, i]) .* blue.(j[:, i]) + ad_blue[:, i] .* (blue.(j[:, i+1]) - blue.(j[:, i]))
    end

    return j

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

using Images, FileIO
img = load("C:\\Users\\Afonso\\Documents\\GitHub\\final_project_fpi\\statue.png")
