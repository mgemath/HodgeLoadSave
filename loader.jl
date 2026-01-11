####################  FAST HODGE DATABASE API  ####################
#
#  Design goals:
#   • Maximum loading performance
#   • One file per (g,n)
#   • Raw integer storage on disk
#   • Single-pass Dict construction
#   • No merge!, no rehash storms
#   • Exact arithmetic via Oscar/Nemo
#
##################################################################

using Serialization
using Oscar   # provides ZZ and QQFieldElem

export hodge_integral
##################################################################
# Types
##################################################################

"""
A key for a Hodge integral:
  (psi-powers, lambda-powers)
"""
const HodgeKey = Tuple{Tuple{Vararg{Int}},Tuple{Vararg{Int}}}

"""
Raw on-disk representation of one Hodge integral.
Only integers are stored.
"""
struct HodgeRaw
  psi::Tuple{Vararg{Int}}
  lambda::Tuple{Vararg{Int}}
  num::Int
  den::Int
end

##################################################################
# Text parsing (only used for conversion)
##################################################################

##################################################################
# Fast loading: ALL (g,n) into one Dict
##################################################################

"""
Load all Hodge integrals with
  g in 1:gMax and n in 1:nMax
from binary .hdg files into a single Dict.
"""
function load_hdg_database(dir::AbstractString,
  gMax::Int,
  nMax::Int)::Dict{HodgeKey,QQFieldElem}

  # --------------------------------------------------------------
  # First pass: load raw blocks and count total size
  # --------------------------------------------------------------
  raw_blocks = Vector{Vector{HodgeRaw}}()
  total_len = 0

  for g in 1:gMax, n in 1:nMax
    file = joinpath(dir, "Hodge_g_$(g) n_$(n).hdg")
    isfile(file) || error("File $file does not exist")

    raws = open(file, "r") do io
      deserialize(io)::Vector{HodgeRaw}
    end

    push!(raw_blocks, raws)
    total_len += length(raws)
  end

  # --------------------------------------------------------------
  # Second pass: build the Dict in one go
  # --------------------------------------------------------------
  H = Dict{HodgeKey,QQFieldElem}()
  sizehint!(H, total_len)

  for raws in raw_blocks
    for r in raws
      H[(r.psi, r.lambda)] = ZZ(r.num) // ZZ(r.den)
    end
  end

  return H
end

function hodge_integral(g::Int64, n::Int64, psi::Vector{Int64}, lambda::Vector{Int64}, H::Dict{HodgeKey, QQFieldElem})::QQFieldElem
  @req length(psi) == n "psi must have length n"
  @req length(lambda) == g "lambda must have length g"
  if sum(psi) + sum((1:g) .* lambda) != 3*g - 3 + n
    return zero(QQ)
  end
  return H[(sort(psi; rev=true)...,), (lambda...,)]
  # return get(H, ((sort(psi; rev=true)...,), (lambda...,)), QQ(0))
end