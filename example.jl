##################################################################
# Example usage (comment out if used as a module)
##################################################################
using BenchmarkTools
include(joinpath(@__DIR__, "converter.jl"))
dir = "$(pwd())/Database/converted"

H = load_hdg_database(dir, 2, 6)
# val = H[((0,1,2), (1,1))]
# 350.828 μs (1690 allocations: 72.67 KiB)
##################################################################
H
hodge_integral(2, 6, [2,1,1,1,1,0], [1,1], H)
# 350.828 μs (1690 allocations: 72.67 KiB)
for i in 1:6
  text_file = pwd()*"/Database/Hodge_g_2 n_$(i).txt"
  hdg_file = pwd()*"/Database/converted/Hodge_g_2 n_$(i).hdg"
  convert_text_to_hdg(text_file, hdg_file)
end
