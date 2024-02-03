#Wrap everything with a module to allow redefinition of type
module RTP
    include("lattice_functions.jl")
    dirs = [[0,1],[0,-1],[1,0],[-1,0]]

    #parameters of the model
    mutable struct Param
        alpha::Float64  #flipping rate
        v::Float64      # hopping rate
        W::Float64      #an arbitrary upper bound of move rate
        Lx::Int64       #system size along x
        Ly::Int64       #system size along y
        N::Int64        #number of particles
        rho_0::Float64  #number of particles per site in average
        rho_m::Float64  #scale of number of particles per site
    end

    #constructor
    function Param(v, alpha, Lx, Ly, rho_0, rho_m)
        W = v + alpha # max possible rate

        N = round(Lx*Ly*rho_0)     #the number of particles

        param = Param(alpha, v, W, Lx, Ly, N, rho_0, rho_m)
        return param
    end

    #arrays containing the position and orientation of the particles, density, and magnetization fields
    mutable struct State
        t::Float64              #time
        pos::Array{Int64,2}     #position of the particle, [particle index, dimension]
        s::Array{Int8, 2}       #orientation of the particle, [0,\pm 1] or [\pm 1, 0]
        rho::Array{Int64,2}     #occupation field
    end

    #constructor
    function State(t, pos, s, Lx, Ly)
        rho = zeros(Int64, Lx, Ly)    #density fields
        N = length(pos[:,1])        #the number of particles

        for n in 1:N    #for each particle
            x, y = pos[n,1], pos[n,2]   #read its position
            rho[x,y] += 1               #add 1 to the density
        end

        st = State(t, pos, s, rho)
        return st
    end

    #main simulator
    function RTP_update!(st, param, t_run, rng)
        #Initialize 
        N, alpha, v, rho_m = param.N, param.alpha, param.v, param.rho_m
        W, Lx, Ly = param.W, param.Lx, param.Ly
        n, x0, x, y0, y, move = 0, 0, 0, 0, 0, 0    #particle index, initial positions, updated positions, selected move index 

        #list of moves
        #1: flip    2~5: +x -x, +y -y, 6: rejection

        #weight vector for the moves. rate for motions along the y-axis are fixed. 
        #W will be used as a normalization factor for the move probability
        #The move trial will be rejected if none of the first 5 is chosen
        w = [alpha/4, alpha/4, alpha/4, alpha/4,  0   ,  W]
            #              tumble                 hop    do nothing

        t = st.t                #read the current time
        t_end = t + t_run       #calculate the ending time
        dt = 1/(N*W)    #time lapse per each step

        while t≤t_end
            n = rand(rng, 1:N)  #randomly select a spin
            x0, y0 = st.pos[n,1], st.pos[n,2]   #read its position
            x = mod1(x0 + st.s[n,1], Lx)
            y = mod1(y0 + st.s[n,2], Ly)
            # hop
            w[5] = hopping_rate(st.rho[x,y], rho_m, v)

            move = tower_sampling(w, W, rng)    #select

            if move==5  # hop
                st.pos[n,1] = x
                st.pos[n,2] = y
                st.rho[x0,y0] -= 1
                st.rho[x,y] += 1

            elseif move<5  # tumble
                st.s[n,1] = dirs[move][1]
                st.s[n,2] = dirs[move][2]
            end

            ####
            ####

            #update the current time
            t += dt
        end
        #recording the time
        st.t = t
    end
    
    function hopping_rate(rho, rho_m, v)
        # rho of individual species?
        if rho<rho_m
            return v*(1-rho/rho_m)
        else
            return v/rho_m
        end
    end
    
end


function make_movie!(st, param, t_gap, n_frame, rng, in_fps)
    #initialize
    println("Movie production started")
    st.t = 0

    prog = 10
    println("Initialization completed, simulation started")
    #animation macro
    anim = @animate for frame in 1:n_frame
        if mod(frame,n_frame÷10)==0   
             println(string(prog)*"% done")
             prog += 10
        end

        RTP.RTP_update!(st, param, t_gap, rng)

        rho_0 = param.rho_0
        rho_m = param.rho_m
        obj = st.rho
        t = st.t
        title = @sprintf("t = %d", t)
        crange = (0,3*rho_m)
        cmap = cgrad(:inferno)

        x_range = range(1, param.Lx, length = param.Lx)
        y_range = range(1, param.Ly, length = param.Ly)

        #heatmap(obj, xlims = (0, param.Lx), ylims = (0, param.Ly), clims = crange, aspect_ratio = 1, label = L"t="*string(st.t[1]))
        heatmap(x_range, y_range, transpose(obj), c = cmap, title=title,
        xlabel = "x", ylabel = "y", xlims = (0, param.Lx), ylims = (0, param.Ly), clims = crange, aspect_ratio = 1)
    end

    println("Simulation complete, producing movie")
    name = "RTP.gif"
    gif(anim, name, fps = in_fps)
end
