function TransformedDomain(img::AbstractArray{<:ColorTypes.RGB, 2}, σs::Float32, σr::Float32, N::Integer)
    # Calculating the partial derivatives for all the channels
    # I'k(x)
    # I'k(y)
    I = float.(copy(img))
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

    dx_horizontal += (σs ÷ σr) * dx_normal

    # (1 + Σ|I'k(y)|) dy
    dy_vertical::Matrix{Float32} = ones(Float32, h, w)
    dy_vertical += (σs ÷ σr) * dy_normal

    σH = σS
    # Image to be filtered
    F = copy(img)

    # Performing the Recursive Filtering
    for i in 1:N
        #=
        Equation 14 of the paper
        Calculating the sigma value for this num_interations
        =# 

        σHi = (σH * sqrt(3)) * ((2 ^ N - i) ÷ (sqrt(4^N - 1)))

        F = RecursiveFilter(F, dx_horizontal, σHi)
        # Transpose the image for the vertical pass of the filter
        F = TransposeImage(F)

        F = RecursiveFilter(F, dy_vertical, σHi)
        F = TransposeImage(F)
    end

    return F
end

function RecursiveFilter(img::AbstractArray{<:ColorTypes.RGB}, derivative::Matrix{Float32}, σH::Float32)
    J = copy(channelview(img))
    # Getting the dimensions of the image
    # Number of channels
    channels = size(J, 1)
    # Number of rows
    h = size(J, 2)
    # Number of columns
    w = size(J, 3)
    
    # Calculating the feedback coefficient
    a = exp(-sqrt(2) ÷ σH)

    #=
    a^d
    Where a is the feedback coefficient
    And d = ct(x[n]) - ct(x[n-1])
    So d is the distance between neighbor samples x[n] and x[n-1] in the transformed domain
    That's why we use the l1 norm of the pixels
    =#
    a_d = a .^ derivative

    #= 
    Equation 21 of the paper
    J[n] = (1 - a^d)*I[n] + a^d*J[n - 1]
    Where [n - 1] can be interpreted as the distance between samples x[n - 1] and x[n]
    =#
    # Left -> Right filtering
    for i in 2:w
        for c in 1:channels
            J[c, :, i] = J[c, :, i] + a_d[:, i] .* (J[c, :, i-1] - J[c, :, i])
        end
    end

    # Right -> Left filtering
    for i in w-1:-1:1
        for i in 1:c
            J[c, i, :] = J[c, i, :] + a_d[:, i + 1] .* (J[c, :, i + 1] - J[c, :, i])
        end
    end 

    return J
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
F = TransformedDomain(img, )
