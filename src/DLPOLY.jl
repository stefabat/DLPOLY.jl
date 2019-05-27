
module DLPOLY

using Printf
using LinearAlgebra

export xyz2dlp,dlp2xyz
export property,rdf,zden
export history
export fragments_dist
export read_CONFIG


"""Convert an `xyz` file into a `CONFIG` file"""
function xyz2dlp(xyzfile,configfile::String="CONFIG")
    # read xyz file line by line, where every line is a string
    xyz = readlines(xyzfile)

    # extract info from first two lines from xyz file
    natms  = parse(Int,xyz[1])
    title = xyz[2]

    # DL POLY variables
    levcfg = 0
    imcon  = 1
    cell   = [25.0  0.0  0.0;
               0.0 25.0  0.0;
               0.0  0.0 25.0]

    # open output file
    fout = open(configfile,"w")

    # print header
    @printf(fout,"%-80s\n",title)
    # print keywords
    @printf(fout,"%10i%10i%10i\n",levcfg,imcon,natms)
    # print simulation cell
    @printf(fout,"%20.10f%20.10f%20.10f\n",cell[1,1],cell[1,2],cell[1,3])
    @printf(fout,"%20.10f%20.10f%20.10f\n",cell[2,1],cell[2,2],cell[2,3])
    @printf(fout,"%20.10f%20.10f%20.10f\n",cell[3,1],cell[3,2],cell[3,3])

    # print atoms
    for i = 1:natms
        line = split(xyz[i+2])
        @printf(fout,"%-8s%8i\n",line[1],i)
        coords = map(x->parse(Float64,x),line[2:4])
        @printf(fout,"%20.10f%20.10f%20.10f\n",coords[1],coords[2],coords[3])
    end
    close(fout)
end


"""Convert a `REVCON` or `CONFIG` file to an `xyz` file"""
function dlp2xyz(CONFIG,xyzfile::String="CONFIG.xyz")
    fin   = open(CONFIG,"r")
    fout  = open(xyzfile,"w")
    title = readline(fin)
    info  = split(readline(fin))
    jump  = parse(Int,info[1])
    bc    = parse(Int,info[2])
    if bc > 0
        box = [map(x->parse(Float64,x),split(readline(fin))) for i=1:3]
    end
    system = readlines(fin)
    close(fin)
    natoms = Int((length(system)/(2+jump)))
    write(fout,"$natoms\n")
    write(fout,title,"\n")
    for i = 1:2+jump:length(system)
        atom   = split(system[i])[1]
        coords = map(x->parse(Float64,x),(split(system[i+1])))
        s = @printf(fout,"%-4s%12.6f%14.6f%14.6f\n",atom,coords[1],coords[2],coords[3])
    end
    close(fout)
end



function read_CONFIG(CONFIG)
    fin   = open(CONFIG,"r")
    title = readline(fin)
    info  = split(readline(fin))
    jump  = parse(Int,info[1])
    bc    = parse(Int,info[2])
    if bc > 0
        box = [map(x->parse(Float64,x),split(readline(fin))) for i=1:3]
    end
    system = readlines(fin)
    close(fin)
    natoms = Int((length(system)/(2+jump)))
    coords = zeros(natoms,3)
    atoms = Vector{String}(undef,natoms)
    j = 1
    for i = 1:2+jump:length(system)
        atom   = split(system[i])[1]
        xyz = map(x->parse(Float64,x),(split(system[i+1])))
        coords[j,1] = xyz[1]
        coords[j,2] = xyz[2]
        coords[j,3] = xyz[3]
        atoms[j] = atom
        j += 1
    end 

    return atoms,coords
end



struct Frame
    number::Int
    atoms::Vector{String}
    coords::Matrix{Float64}
end



"""Read a `HISTORY` file"""
function history(HISTORY)
    fin = open(HISTORY,"r")
    # skip header line
    readline(fin)
    # read info on second line
    info = split(readline(fin))
    keytrj  = parse(Int,string(info[1]))
    natms   = parse(Int,string(info[3]))
    nframes = parse(Int,string(info[4]))
    if keytrj > 0
        error("possible to read a `HISTORY` file with positions only")
    end

    # initialize some arrays
    step   = Vector{Int}(undef,nframes)
    time   = Vector{Float64}(undef,nframes)
    frames = Vector{Frame}(undef,nframes)


    # start loop over frames
    for i = 1:nframes
        # read record i
        reci = split(readline(fin))
        step[i] = parse(Int,string(reci[2]))
        time[i] = parse(Float64,string(reci[7]))
        
        # skip 3 lines, i.e., records ii, iii, iv
        readline(fin);readline(fin);readline(fin)
        
        # now we read all atoms and the coords to
        # create the Frame and save it in the frames
        # array, but before I initialize some shitty
        # temporay arrays because I don't want to
        # think how to do it more elegantly
        atm = Vector{String}(undef,natms)
        xyz = Matrix{Float64}(undef,(natms,3))
        
        for j = 1:natms
            # read record a
            reca = split(readline(fin))
            # save atom label
            atm[j] = reca[1]

            # read record b, i.e. positions
            recb = map(string,split(readline(fin)))
            # save coordinates of that atom
            xyz[j,:] = map(x->parse(Float64,x),recb)
        end

        # create Frame
        frames[i] = Frame(i,atm,xyz)

    end

    return frames
end


function fragments_dist(mol1::Matrix, mol2::Matrix)
    natm1 = size(mol1,1)
    natm2 = size(mol2,1)

    min_dist = maxintfloat()*ones(natm1)

    for i=1:natm1
        for j=1:natm2
            d = norm(mol1[i,:]-mol2[j,:])
            if d <= min_dist[i]
                min_dist[i] = d
            end
        end
    end

    return min_dist
end


"""Extract the radial distribution function from an `RDFDAT` file"""
function rdf(RDFDAT)
    fin = open(RDFDAT,"r")
    title = readline(fin)
    info = map(x->parse(Int,x),split(readline(fin)))
    rdf = Matrix(info[2]+1,info[1]+1)
    for j = 1:info[1]
        for i = 1:info[2]+1
            if j > 1 && i > 1
                rdf[i,j+1] = parse(Float64,split(readline(fin))[2])
            elseif j == 1 && i > 1
                rdf[i,1:2] = map(x->parse(Float64,x),split(readline(fin)))
            elseif j > 1 && i == 1
                label = split(readline(fin))
                rdf[i,j+1] = string(label[1],"-",label[2])
            elseif j == 1 && i == 1
                rdf[i,j] = "R"
                label = split(readline(fin))
                rdf[i,j+1] = string(label[1],"-",label[2])
            end
        end
    end
    return rdf
end


"""Extract the Z density profile from an `ZDNDAT` file"""
function zden(ZDNDAT)
    fin = open(ZDNDAT,"r")
    title = readline(fin)
    info = map(x->parse(Int,x),split(readline(fin)))
    zden = Matrix(info[2]+1,info[1]+1)
    for j = 1:info[1]
        for i = 1:info[2]+1
            if j > 1 && i > 1
                zden[i,j+1] = parse(Float64,split(readline(fin))[2])
            elseif j == 1 && i > 1
                zden[i,1:2] = map(x->parse(Float64,x),split(readline(fin)))
            elseif j > 1 && i == 1
                label = split(readline(fin))
                zden[i,j+1] = string(label[1])
            elseif j == 1 && i == 1
                zden[i,j] = "Z"
                label = split(readline(fin))
                zden[i,j+1] = string(label[1])
            end
        end
    end
    return zden
end


"""Extract `record` from `STATIS` file"""
function property(STATIS, record::String)
    f = open(STATIS,"r")
    title = readline(f)
    units = split(readline(f),"=")[2]
    # read the rest of the file entirely
    data = readlines(f)
    close(f)

    i = 1
    step = Int[]
    prop = Float64[]

    # check which property is asked
    if record == "engcns" || record == "energy"
        loffset = 1     # line of the block
        lentry  = 1     # position on the line
    elseif record == "temp" || record == "temperature"
        loffset = 1     # line of the block
        lentry  = 2     # position on the line
    elseif record == "engcfg"
        loffset = 1     # line of the block
        lentry  = 3     # position on the line
    elseif record == "engsrp"
        loffset = 1     # line of the block
        lentry  = 4     # position on the line
    elseif record == "engcpe"
        loffset = 1     # line of the block
        lentry  = 5     # position on the line
    elseif record == "rottemp"
        loffset = 3     # line of the block
        lentry  = 1     # position on the line
    elseif record == "vir" || record == "virial"
        loffset = 3     # line of the block
        lentry  = 2     # position on the line
    elseif record == "virsrp"
        loffset = 3     # line of the block
        lentry  = 3     # position on the line
    elseif record == "vircpe"
        loffset = 3     # line of the block
        lentry  = 4     # position on the line
    elseif record == "volume"
        loffset = 4     # line of the block
        lentry  = 4     # position on the line
    elseif record == "pressure"
        loffset = 6     # line of the block
        lentry  = 2     # position on the line
    else
        error("not implemented yet")
    end

    # loop over blocks
    while i < length(data)
        # info[1] = step, info[2] = time, info[3] = number of entries per block
        #info = map(x->parse(Int,x),split(data[i]))
        info = Vector(undef,3)
        info[1] = parse(Int,split(data[i])[1])
        info[2] = parse(Float64,split(data[i])[2])
        info[3] = parse(Int,split(data[i])[3])
        # save step number
        push!(step,info[1])
        # determine the blocksize
        bs = Int(ceil(info[3]/5))
        # temperature is always in the first row of the block
        line = map(x->parse(Float64,x),split(data[i+loffset]))
        # second entry corresponds to temperature
        push!(prop,line[lentry])
        # update iterator to jump to the next block
        i += bs + 1
    end

    return step,prop
end

end
