function get_correlation!(map, param, cor_map)
    #get density density and current current correlation

    prod_conv = zeros(ComplexF64, param.Lx, param.Ly)

    #density correlation
    mapk = fft(map)
    for x in 1:param.Lx, y in 1:param.Ly
        prod_conv[x,y] = mapk[x,y]*mapk[mod1(2-x, param.Lx), mod1(2-y,param.Ly)]
    end

    copyto!(cor_map, real.(ifft(prod_conv)/(param.Lx*param.Ly)) )
end

#The function for randomly selecting an option according to the weight. w_sum is the normalization factor
function tower_sampling(weights, w_sum, rng)
    key = w_sum*rand(rng)

    selected = 1
    gathered = weights[selected]
    while gathered < key
        selected += 1
        gathered += weights[selected]
    end

    return selected
end

#Calculate the histogram array of a 2d array
function density_histogram(rho)
    state_flat = vec(rho)
    max_rho = maximum(state_flat)
    rho_histogram = zeros(Int64, max_rho+1)

    for i in eachindex(state_flat)
        rho_histogram[state_flat[i]+1] += 1
    end

    return rho_histogram
end