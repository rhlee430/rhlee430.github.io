---
title: Exact Diagonalization of Spin Models
---
## Indexing bases

Basis states of spin-$l$ models are represented with a string of 0 and 1. We should map this string to an integer by

$$
 |s_1, s_2, \dots, s_N\rangle \rightarrow |a\rangle \;
 \text{with} \; a={\sum_{i=1}^N s_i \times (2l+1)^i}.
$$
$l=1/2$ case appears in most cases. For those cases, labeling a basis with an integer $a$ is equivalent to regarding the string as a binary number. Thus, we could make functions that transform an integer(or a string) index to the other side.

```jl
 ## mapping between integer index and bit string
function idx_to_bitString(N::Int64, idx::Int64)::Vector{Bool}

	bitString = zeros(Bool, N)
	
	for i in 0:N-1
		bitString[1+i] = idx & 1
		idx >>= 1
	end

	return bitString
end

function bitString_to_idx(N::Int64, bitString::Vector{Bool})::Int64

	idx = 0

	for i in 0:N-1
		if bitString[1+i]
			idx += 1<<i
		end
	end

	return idx
end
```

We will now consider 1D transverse field Ising model to show how to construct the Hamiltonian in a matrix form.
## 1D Transverse Field Ising Model (TFIM)

The model Hamiltonian is 
$$
H = -{\sum_{\langle ij\rangle} X_iX_j} \; - h{\sum_{i}Z_i}.
$$
As we now how $H$ acts on $|\{s_i\}\rangle$, we can interpret it with $|a\rangle$ and fill out components of the Hamiltonian matrix.

```jl
# Example: 1D TFIM
# construct dense matrix representation of many-body Hamiltonian
function construct_matH_dense(siteNum::Int64, h::Float64; pbc::Bool=true)::Hermitian{Matrix{Float64}}

	# combined Hilbert space dimension
	matDim = 1 << siteNum

	# construct matrix
	matH = zeros(Float64, matDim, matDim)

	for idx in 0:matDim - 1
		bitstring = idx_to_bitString(siteNum, idx)
		
		# XX term
		for i in 0:siteNum - 1
		
			# (right-)nearest neighbor site
			if !pbc && i==siteNum-1; continue; end
			ip = (i+1)%siteNum
			
			# X[i]X[ip]|idx> = |idx_>
			bitString_ = copy(bitString)
			bitString_[1+i] = !bitString_[1+i]
			bitString_[1+ip] = !bitString_[1+ip]
			idx_ = bitString_to_idx(siteNum, bitString_)
			matH[1+idx_, 1+idx] += -1.0
		end
		
		# Z term
		for i in 0:siteNum - 1
			
			# Z[i]|idx> = (val)*|idx>
			# where (val) = +h or -h
			if bitString[1+i] == false
				matH[1+idx, 1+idx] += -h
			else
				matH[1+idx, 1+idx] += +h
			end
		end
	end

	return Hermitian(matH)
end
```

For the Hamiltonian with complex components, we should use *Complex64* instead of *Float64*.