using JLD2
using MAT

function move_to_mat(mom, tra; prefix = "averages")
    file = jldopen("channel_$(prefix)_$(mom * 10 + tra).jld2", "r")

    iter = keys(file["timeseries/t"])[end]
    Pu = file["timeseries/Pv/$iter"]
    Pv = file["timeseries/Pu/$iter"]
    Pw = file["timeseries/Pu/$iter"]
    Bx = file["timeseries/Pu/$iter"]
    By = file["timeseries/Pu/$iter"]
    Bz = file["timeseries/Pu/$iter"]

    U = file["timeseries/Uⁿ⁻¹/$iter"]
    V = file["timeseries/Vⁿ⁻¹/$iter"]
    W = file["timeseries/Wⁿ⁻¹/$iter"]
    B = file["timeseries/bⁿ⁻¹/$iter"]

    matfile = matopen("channel_$(prefix)_$(mom * 10 + tra).mat", "w")
    write(matfile, "Pu", Pu)
    write(matfile, "Pv", Pv)
    write(matfile, "Pw", Pw)
    write(matfile, "Bx", Bx)
    write(matfile, "By", By)
    write(matfile, "Bz", Bz)
    write(matfile, "U", U)
    write(matfile, "V", V)
    write(matfile, "W", W)
    write(matfile, "B", B)
    close(matfile)
    
    return nothing
end