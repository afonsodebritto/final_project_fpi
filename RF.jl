function TransformedDomainRecursiveFilter(img, σs, σr, N)
    # Calculating the partial derivatives for all the channels
    # I'k(x)
    # I'k(y)
    I = float.(channelview(img))
    J = I
    dx_partial = diff(J, dims=3)
    dy_partial = diff(J, dims=2)

    # Creating the accumulated derivative matrixes for all the channels
    channels = size(J, 1)
    h = size(J, 2)
    w = size(J, 3)
    dx_normal = zeros(h, w)
    dy_normal = zeros(h, w)
    # Calculating the l1-normal distance of the pixels
    # Σ|I'k(x)|
    # Where k is the color channel
    for c in 1:channels
        for i in 1:h
            for j in 2:w
                dx_normal[i, j] = dx_normal[i, j] + abs(dx_partial[c, i, j - 1])
            end
        end
    end

    # Σ|I'k(Y)|
    # Where k is the color channel
    for c in 1:channels
        for i in 2:h
            for j in 1:w
                dy_normal[i, j] = dy_normal[i, j] + abs(dy_partial[c, i - 1, j])
            end
        end
    end

    # Calculating the horizontal and vertical derivatives of the transformed domain
    # (1 + (σs÷σr)*Σ|I'k(x)|) dx    
    dx_horizontal = ones(h, w)
    for i in 1:h
        for j in 1:w
            dx_horizontal[i, j] = dx_horizontal[i, j] + (σs ÷ σr) * dx_normal[i, j]
        end
    end

    # (1 + Σ|I'k(y)|) dy
    dy_vertical = ones(h, w)
    for i in 1:h
        for j in 1:w
            dy_vertical[i, j] = dy_vertical[i, j] + (σs ÷ σr) * dy_normal[i, j]
        end
    end 

    σH = σs
    # Image to be filtered
    F = I

    # The vertical pass is performed on a transposed image, therefore we transpose the vertical derivative
    dy_vertical = TransposeMatrix(dy_vertical)

    # Performing the Recursive Filtering
    for i in 1:N
        #=
        Equation 14 of the paper
        Calculating the sigma value for this num_interations
        =# 

        σHi = σH * sqrt(3) * 2^(N - i) / sqrt(4^N - 1)

        F = RecursiveFilter(F, dx_horizontal, σHi)
        # Transpose the image for the vertical pass of the filter
        F = TransposeImage(F)

        F = RecursiveFilter(F, dy_vertical, σHi)
        # Transpose the image back to the original orientation
        F = TransposeImage(F)
    end

    return F
end

function RecursiveFilter(img, derivative, σH)
    J = img
    # Getting the dimensions of the image
    # Number of channels
    channels = size(J, 1)
    # Number of rows
    h = size(J, 2)
    # Number of columns
    w = size(J, 3)
    
    # Calculating the feedback coefficient
    a = exp(-sqrt(2) / σH)
    a_d = fill(a, (h, w))

    #=
    a^d
    Where a is the feedback coefficient
    And d = ct(x[n]) - ct(x[n-1])
    So d is the distance between neighbor samples x[n] and x[n-1] in the transformed domain
    That's why we use the l1 norm of the pixels
    =#

    for i in 1:h
        for j in 1:w
            a_d[i, j] = a ^ derivative[i, j]
        end
    end

    #= 
    Equation 21 of the paper
    J[n] = (1 - a^d)*I[n] + a^d*J[n - 1]
    Where [n - 1] can be interpreted as the distance between samples x[n - 1] and x[n]
    =#
    # Left -> Right filtering
    for c in 1:channels
        for i in 1:h
            for j in 2:w
                J[c, i, j] = clamp(J[c, i, j] + a_d[i, j] * (J[c, i, j-1] - J[c, i, j]), 0, 1)
            end
        end
    end

    # Right -> Left filtering
    for i in 1:h
        for j in w-1:-1:1
            for c in 1:channels
                J[c, i, j] = clamp.((J[c, i, j]) + a_d[i, j+1] * (J[c, i, j+1] - J[c, i ,j]), 0, 1)
            end
        end
    end

    return J
end

function TransposeMatrix(matrix)
    h = size(matrix, 1)
    w = size(matrix, 2)

    # Creating a transposed matrix with the correct size
    transposed_matrix = similar(matrix, w, h)

    # Assigning each row of the original matrix to a column in the transposed one
    for i in 1:h
        transposed_matrix[:, i] = matrix[i, :]
    end

    return transposed_matrix
end

function TransposeImage(img)
    # Getting the height of the image
    h = size(img, 2)
    
    # Creating a transposed image with the correct size
    transposed_image = similar(img, size(img, 1), size(img, 3), size(img, 2))
    
    # Assigning each row of the original img to a column in the transposed one
    for i in 1:h
        transposed_image[:, :, i] = img[:, i, :]
    end
    
    return transposed_image
end


using Images, FileIO
img_original = load("C:\\Users\\Afonso\\Documents\\GitHub\\final_project_fpi\\statue.png")
img_filtrada = colorview(RGB, TransformedDomainRecursiveFilter(img_original, 60, 0.4, 3))
mosaic(img_original, img_filtrada; nrow = 1)