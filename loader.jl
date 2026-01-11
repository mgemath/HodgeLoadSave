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
# const HodgeKey = Tuple{Vector{Int},Vector{Int}}

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

@inline function hodge_integral(
  g::Int,
  n::Int,
  C::Vector{Int},
  valG::Int,
  Ev::Int,
  markPsis::Vector{Int},
  lambda::Vector{Int},
  H::Dict{HodgeKey,QQFieldElem},
  psi_buf::Vector{Int}   # preallocated buffer psi_buf = Vector{Int}(undef, n)   # n = Ev + length(markPsis)
)::QQFieldElem

  # @boundscheck begin
  #   length(lambda) == g || error("lambda length mismatch")
  #   length(psi_buf) == n || error("psi buffer length mismatch")
  # end

  # build psi without allocations
  @inbounds begin
    k = 1
    for i in valG+1:valG+Ev
      psi_buf[k] = C[i]
      k += 1
    end
    for x in markPsis
      psi_buf[k] = x
      k += 1
    end
  end

  # dimension check (cheap)
  s = 0
  @inbounds for x in psi_buf
    s += x
  end
  @inbounds for i in 1:g
    s += i * lambda[i]
  end
  s == 3g - 3 + n || return zero(QQ)

  # canonical key
  sort!(psi_buf; rev=true)
  key = (Tuple(psi_buf), Tuple(lambda))

  return get(H, key, zero(QQ))
end


function _check_file(dir::AbstractString,
  gMax::Int,
  nMax::Int)::Bool

  for g in 1:gMax, n in 1:nMax
    file = joinpath(dir, "Hodge_g_$(g) n_$(n).hdg")
    isfile(file) || return false
  end
  return true

end