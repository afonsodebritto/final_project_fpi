function TransformedDomainInterpolatedConvolution(img, σs, σr, N)
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

    # ∫(1 + Σ|I'k(x)|) dx
    ct_H = cumsum(dx_horizontal, dims=2)

    # ∫(1 + Σ|I'k(y)|) dy
    ct_V = cumsum(dy_vertical, dims=1)

    # The vertical pass is performed on a transposed image, therefore we transpose the vertical derivative
    ct_V = TransposeMatrix(ct_V)

    # Performing the Normalized Convolution
    for i in 1:N
        #=
        Equation 14 of the paper
        Calculating the sigma value for this num_interations
        =# 

        σHi = σH * sqrt(3) * 2^(N - i) / sqrt(4^N - 1)
        box_filter_radius = σHi * sqrt(3)

        F = NormalizedConvolution(F, ct_H, box_filter_radius)
        # Transpose the image for the vertical pass of the filter
        F = TransposeImage(F)

        F = NormalizedConvolution(F, ct_V, box_filter_radius)
        # Transpose the image back to the original orientation
        F = TransposeImage(F)
    end

    return F
end

function InterpolatedConvolution(img, domain_transform, box_filter_radius)
    h = size(img, 2)
    w = size(img, 3)
    channels = size(img, 1)

    # Calculating the limits
    lower_limit = domain_transform .- box_filter_radius
    upper_limit = domain_transform .+ box_filter_radius

    # Calculating the indexes of the limit pixels
    lower_indexes = zeros(Int, size(domain_transform))
    upper_indexes = zeros(Int, size(domain_transform))
    
    for i in 1:h
        domain_transform_cur_row = domain_transform[i, :]
        push!(domain_transform_cur_row, ∞)
    
        lower_limit_cur_row = lower_limit[i, :]
        upper_limit_cur_row = upper_limit[i, :]
    
        lower_limit_cur_row_index = zeros(Int, 1, w)
        upper_limit_cur_row_index = zeros(Int, 1, w)
        
        # Finding the lower and upper indexes for the first element
        lower_limit_cur_row_index[1] = findfirst(domain_transform_cur_row .> lower_limit_cur_row[1])
        upper_limit_cur_row_index[1] = findfirst(domain_transform_cur_row .> upper_limit_cur_row[1])
    
        for j = 2:w
            lower_limit_cur_row_index[j] = lower_limit_cur_row_index[j-1] + 
                findfirst(domain_transform_cur_row[lower_limit_cur_row_index[j-1]:end] .> lower_limit_cur_row[j]) - 1
            
            upper_limit_cur_row_index[j] = upper_limit_cur_row_index[j-1] + 
                findfirst(domain_transform_cur_row[upper_limit_cur_row_index[j-1]:end] .> upper_limit_cur_row[j])  - 1
        end

        lower_indexes[i,:] = lower_limit_cur_row_index
        upper_indexes[i,:] = upper_limit_cur_row_index
    end

    # Calculating the heights of the trapezoids
    heights = 0.5 .* (img[:,:,2:end] + img[:,:,1:end-1])
    # Calculating the width of the trapezoids
    widths = domain_transform[:,2:end] - domain_transform[:,1:end-1]
    # Calculating the areas of the trapezoids
    areas = similar(heights)
    for i in 1:h
        for j in 1:w-1
            for c in 1:channels
                areas[c,i,j] = heights[c,i,j] * widths[i,j]
            end
        end
    end 

    # Calculating the summed area table
    summed_area_table = zeros(channels, h, w)
    summed_area_table[:,:,2:end] = cumsum(areas, dims=3)

    # Adding paddings
    # Padding the image
    padded_img = zeros(channels,h,w+2)
    padded_img[:,:,2:w+1] = img[:,:,:]
    padded_img[:,:,1] = img[:,:,1]
    padded_img[:,:,w+2] = img[:,:,w]
    # Padding the summed area table
    padded_summed_area_table = zeros(channels, h, w+2)
    padded_summed_area_table[:,:,2:w+1] = summed_area_table[:,:,:]
    # Padding the domain transform coordinates
    padded_domain_transform = zeros(h,w+2)
    padded_domain_transform[:,2:w+1] = domain_transform[:,:]
    padded_domain_transform[:,1] = domain_transform[:,1]
    padded_domain_transform[:,w+2] = domain_transform[:,w]

    # Assuming pixels out of the bounds are equal to nearest border value
    for i in 1:w+2
        padded_domain_transform[i,1] = padded_domain_transform[i,1] - 1.2 * box_filter_radius
    end
    for i in 1:w+2
        padded_domain_transform[i,w+2] = padded_domain_transform[i,w+2] - 1.2 * box_filter_radius
    end

    lower_indexes = lower_indexes .+ 1

    F = zeros(size(img))


    # Filtering the image
    for i in 1:h
        for j in 1:w
            for c in 1:channels
                F[c, i, j] = (summed_area_table[c, i, upper_indexes[i,j]] - summed_area_table[c, i, lower_indexes[i,j]]) / (upper_indexes[i,j] - lower_indexes[i,j])
            end
        end
    end

    return F

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

using Images, FileIO, Infinity, ImageFiltering
img = load("C:\\Users\\Afonso\\Documents\\GitHub\\final_project_fpi\\statue.png")
img_filtrada = colorview(RGB, TransformedDomainNormalizedConvolution(img, 60, 0.4, 3))