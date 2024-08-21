include("src/visualize.jl")

const Nx = 200
const Ny = 400
const Nz = 90

function read_variable(file, Nvar)
     var = Array{Float32}(undef, Nx*Ny*Nz*Nvar)
     read!(file, var)
     var = bswap.(var) |> Array{Float64}
     var = reshape(var, Nx, Ny, Nz, Nvar)
     var = reverse(var, dims = 3)

     return var
end

