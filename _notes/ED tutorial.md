---
title: Exact Diagonalization
---
In exact diagonalization (ED), we construct many-body Hamiltonian to a matrix form, and by diagonalizing it, we obtain full spectrum of the Hamiltonian. In this note, we deal with ED calculations with *Julia*. 
## Interacting Spins

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

$$H = -{\sum_{\langle ij\rangle} X_iX_j} - h{\sum_{i}Z_i}.$$

As we now how $H$ acts on spins, we can interpret it with other labels and fill out components of the Hamiltonian matrix.

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

For the Hamiltonian with complex components, we should use *ComplexF64* instead of *Float64*.

## Fermions

An indistinguishability of fermions make numerics a bit complicated. Due to second quantization, we could denote bases of fermionic Fock space with occupation numbers of each sites (including spin d.o.f.). By the Pauli exclusion principle, 0 or 1 is only values for occupations. Moreover, if we define an occupation number of the generic site index $i$ as $n_i$, we can make a correspondence via

$$
|\{n_i\}\rangle \Rightarrow a = \sum_{i=0}^{N-1}n_i 2^i \, \text{th vector},
$$

where $N$ is a dimension of single particle Hilbert space, i.e., (# of sites) $\times$ (spin d.o.f.).

We will now apply this formalism to the one-dimensional Hubbard model.

## 1D Hubbard Model

The Hamiltonian of the model is

$$
H = -t\sum_{\langle i j \rangle, \sigma} \left( c_{i, \sigma}^\dagger c_{j, \sigma} + h.c. \right) + U \sum_i n_{i\uparrow} n_{i\downarrow} - \mu \sum_{i, \sigma} n_{i, \sigma}.
$$

Following constructs the electron Hubbard model Hamiltonian in a matrix form.
```jl
# construct dense matrix representation of the Hubbard Hamiltonian
function construct_matH_dense(parms::Tuple{Float64, Float64, Float64}, siteNum::Int64 ; pbc::Bool= false)::Hermitian{ComplexF64, Matrix{ComplexF64}}

	# parameters
	t, U, mu = parms
	
	# combined Hilbert space dimension
	totNum = 2 * siteNum
	matDim = 1 << totNum
	
	# construct matrix
	matH = zeros(ComplexF64, matDim, matDim)
	for idx in 0:matDim - 1
		sqState = idx_to_bitString(totNum, idx)
	  
		# hopping term
		for i in 0:siteNum - 1, s in 0:1
			
			# (right-)nearest neighbor siteNum
			if !pbc && i == siteNum - 1
				continue
			end
			ip = (i+1) % siteNum
				
			# c+[ip, s] c[i, s] |idx> = (val)*|idx_>
			# or c[ip, s] c+[i, s] |idx> = (val)*|idx_>
			if sqState[2*i+s+1] == sqState[2*ip+s+1]
				continue
			end
			
			sign_JW = 1.0
			if sum(@view sqState[1:2*i+s])%2 == 1; sign_JW = - sign_JW; end
			sqState[2*i+s+1] = !sqState[2*i+s+1]
			
			if sum(@view sqState[1:2*ip+s])%2 == 1; sign_JW = - sign_JW; end
			sqState[2*ip+s+1] = !sqState[2*ip+s+1]
			
			idx_ = bitString_to_idx(totNum, sqState)
			
			matH[1+idx_, 1+idx] += sign_JW * (-t)
			
			# restore to original sqState
			sqState[2*i+s+1] = !sqState[2*i+s+1]
			sqState[2*ip+s+1] = !sqState[2*ip+s+1]
		end
		
		# Hubbard interaction
		for i in 0:siteNum - 1
			if sqState[2*i+1] && sqState[2*i+2]
				matH[1+idx, 1+idx] += U
			end
		end
		
		# chemical potential
		matH[1+idx, 1+idx] += -mu * sum(sqState)
	end
	return Hermitian(matH)
end
```

For example, we could obtain energy eigenvalues of the model with two sites as follows.

```jl
# 2 sites, t=1.0, U=50.0, mu=0.0
matH_dense = construct_matH_dense((1.0, 50.0, 0.0), 2)


# dense diagonalization
vals, vecs = eigen(matH_dense)
println(vals)
```

Low energy energy spectrum of Hubbard model with a half-filling could be explained with the *superexchange* interaction, which favors a spin singlet state.
