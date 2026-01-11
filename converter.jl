include(joinpath(@__DIR__, "loader.jl"))


const LINE_RE = r"""
^
\s*(\d+)\s*,\s*(\d+)\s*,
\s*\(([^)]*)\)\s*,\s*
\s*\(([^)]*)\)\s*;\s*
(\d+)\s*//\s*(\d+)
\s*$
"""x

"""
Parse one line of the text database into a HodgeRaw entry.
"""
function parse_line_raw(line::String)::HodgeRaw
  m = match(LINE_RE, line)
  m === nothing && error("Invalid line format:\n$line")

  g = parse(Int, m.captures[1])
  n = parse(Int, m.captures[2])

  psi_data =
    isempty(strip(m.captures[3])) ?
    Int[] :
    parse.(Int, split(m.captures[3], ',', keepempty=false))

  lambda_data =
    isempty(strip(m.captures[4])) ?
    Int[] :
    parse.(Int, split(m.captures[4], ',', keepempty=false))

  @assert length(psi_data) == n
  @assert length(lambda_data) == g

  HodgeRaw(
    ntuple(i -> psi_data[i], n),
    ntuple(i -> lambda_data[i], g),
    parse(Int, m.captures[5]),
    parse(Int, m.captures[6])
  )
end

##################################################################
# Conversion: text file → binary (.hdg)
##################################################################

"""
Convert a single text file (fixed g,n) into a binary .hdg file.
"""
function convert_text_to_hdg(text_file::AbstractString,
  hdg_file::AbstractString)
  raws = HodgeRaw[]

  open(text_file, "r") do io
    for line in eachline(io)
      _HodgeRaw = parse_line_raw(line)
      # skip zeroes if desired
      # iszero(_HodgeRaw.num) && continue
      push!(raws, _HodgeRaw)
    end
  end

  open(hdg_file, "w") do io
    serialize(io, raws)
  end
end

