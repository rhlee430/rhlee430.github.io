---
title: Exact Diagonalization of Composite Fermi Liquids
---
### Code syntax highlighting

You can add code blocks with full syntax color highlighting by wrapping code snippet in triple backticks and specifying the type of the code (`js`, `rb`, `sh`, etc.):

```jl
function idx_to_bitString(N::Int64, idx::Int64)::Vector{Bool}

bitString = zeros(Bool, N) for i in 0:N-1

bitString[1+i] = idx & 1

idx >>= 1

end

return bitString end

function bitString_to_idx(N::Int64, bitString::Vector{Bool})::Int64

idx = 0  
for i in 0:N-1

if bitString[1+i] idx += 1<<i

end end

return idx end
```