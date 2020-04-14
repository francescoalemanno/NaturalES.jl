import Base.size
import Base.getindex

struct SymProd{T,V,L} <: AbstractMatrix{T}
    vec::V
    function SymProd(x::vec) where {T,vec <: AbstractVector{T}}
        new{T,vec,length(x)}(x)
    end
end

Base.size(x::SymProd{T,V,L}) where {T,V,L} = (L,L)

function Base.getindex(A::SymProd{T,V,L},i,j) where {T,V,L}
    @boundscheck begin
        (1<=i<=L) || throw(BoundsError(A.vec,i))
        (1<=j<=L) || throw(BoundsError(A.vec,j))
    end
    @inbounds A.vec[i]*A.vec[j]
end
