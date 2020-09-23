################################################################################
#################### Scaled Grid Interpolation #################################
################################################################################
# The following code generates three objects:
    # 1) A "TypeScale" type that transforms ranges into a type of scaled ranges
    # 2) A "TypeRange" type that foroms scaled ranges of a given type.
    #   These ranges behave "just as" normal ranges
    # 3) A ScaledInterpolations function that wraps "interpolations.jl" to be used
    #   with TypeRange ranges
# Types of Scales: parameters bounds for grid (a,b) and curvature (θ)
    # p = PolyScale(a,b,θ)
        # Maps a number x∈[0,1] to a polynomial grid with curvature θ: p(x)∈[a,b]
        # Usage: p(x), where x is a number, or p.(x) where x is vector/range
        # p(x) = a + (b - a) * x^θ
    # p = InvPolyScale(a,b,θ)
        # Maps a number y∈[a,b] from a polynomial grid with curvature θ to the unit interval: p(x)∈[0,1]
        # Usage: p(y), where x is a number, or p.(y) where y is vector/range
        # p(y) = ((y-a) / (b-a) )^(1/θ)
    # p = ExpScale(a,b,θ)
        # Maps a number x∈[0,1] to a exponential grid with curvature θ: p(x)∈[a,b]
        # Usage: p(x), where x is a number, or p.(x) where x is vector/range
        # p(x) = a + (b - a) * ( (exp(θ*x) - 1)/(exp(θ) - 1) )
    # p = LogScale(a,b,θ)
        # Maps a number y∈[a,b] from an exponential grid with curvature θ to the unit interval: p(x)∈[0,1]
        # Usage: p(y), where x is a number, or p.(y) where y is vector/range
        # p(y) = log( (y-a)*(exp(θ) - 1)/(b - a) + 1 ) / θ
# Types of Scaled Ranges
    # grid = PolyRange(a,b,θ=t,N=n) = PolyScale(a,b,t)(range(0,1,length=n))
    # grid =  ExpRange(a,b,θ=t,N=n) =  ExpScale(a,b,t)(range(0,1,length=n))
# ScaledInterpolations
    # itp = ScaledInterpolations(grid,f_grid,Itp_Type)
        # grid is either a PolyRange or ExpRange type
        # f_grid is a vector with the values of the function at the grid nodes
        # Itp_Type is a type from the interpolations.jl package
            # Examples: BSpline(Cubic(Line(OnGrid()))) // shortcut: CubicSpline
            #           BSpline(Cubic(Flat(OnGrid())))
            #           FritschButlandMonotonicInterpolation() // shortcut: MonotoneSpline
            #           BSpline(Linear()) // shortcut: LinearInterp
# ScaledInterpolations with multiple dimensions
    # itp_md = ScaledInterpolations( (grid_1,...,grid_N) , (f_grid_1,...,f_grid_N) , (Type_1,...,Type_N)  )     
# Jacob Adenbaum, 2020


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Define Scaled types
    #-------------------------------------------------------------------------------
    # Polynomial Scaling
    abstract type Scale{T}       end
    abstract type ScaledRange{T} end

    for Foo in [:PolyScale, :InvPolyScale]

        @eval struct $Foo{T,TF} <: Scale{T}
            a::T
            b::T
            θ::TF
        end

        @eval function $Foo(a, b, θ = 1)
            T = reduce(promote_type, map(typeof, (a,b)))
            $Foo(convert(T, a), convert(T, b), θ)
        end
    end


    (p::PolyScale)(x)           = p.a + (p.b - p.a) * x^p.θ
    (p::InvPolyScale)(y)        = ((y-p.a) / (p.b-p.a) )^(1/p.θ)
    Base.inv(p::PolyScale)      = InvPolyScale(p.a, p.b, p.θ)
    Base.inv(p::InvPolyScale)   = PolyScale(p.a, p.b, p.θ)

    #-------------------------------------------------------------------------------
    # Exponential Scaling
    for Foo in [:ExpScale, :LogScale]

        @eval struct $Foo{T, TF} <: Scale{T}
            a::T
            b::T
            θ::TF
            et::T
            s::T
        end

        @eval function $Foo(a, b, θ = 1)
            s = (exp(θ) - 1)/(b - a)
            et= exp(θ)
            T = reduce(promote_type, map(typeof, (a,b,θ, et,s)))
            $Foo(convert(T, a), convert(T, b), θ, et,s)
        end
    end

    function (p::ExpScale)(x)
        t  = (exp(p.θ * x) - 1)/(p.et - 1)
        return p.a + (p.b - p.a) * t
    end

    function (p::LogScale)(y)
        θ  = p.θ
        log( (y-p.a)*p.s + 1 ) / p.θ
    end

    Base.inv(p::ExpScale)       = LogScale(p.a, p.b, p.θ)
    Base.inv(p::LogScale)       = ExpScale(p.a, p.b, p.θ)

    #-------------------------------------------------------------------------------
    # Finalize Types
    scales = Dict(:PolyRange => :PolyScale, :ExpRange => :ExpScale)


    Base.inv(::Type{ExpScale}) = LogScale
    Base.inv(::Type{LogScale}) = ExpScale
    Base.inv(::Type{PolyScale})= InvPolyScale
    Base.inv(::Type{InvPolyScale}) = PolyScale

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Define Scaled Ranges types
    #-------------------------------------------------------------------------------
    # Polynomial Range
    for Typ in [:PolyRange] @eval begin
        struct $Typ{T,TF} <: AbstractRange{T}
            a::T
            b::T
            θ::TF
            N::Int
            vals::Vector{Float64}
            scale::$(scales[Typ]){T,TF}
            invscale::inv($(scales[Typ])){T,TF}
        end

        function $Typ(a,b; θ = 1, N = 100)
            T       = reduce(promote_type, map(typeof, (a,b,θ, 1.0)))
            Ta      = convert(T, a)
            Tb      = convert(T, b)

            vals    = domscale($Typ)(Ta,Tb,θ).(LinRange(0, 1, N))
            scale   = $(scales[Typ])(Ta, Tb, θ)
            invscale= inv($(scales[Typ])(Ta, Tb, θ))
            $Typ(convert(T, a), convert(T, b), θ, N, vals, scale, invscale)
        end

        domscale(p::$Typ)      = p.scale
        domscale(::Type{$Typ}) = $(scales[Typ])

        function rangescale(r::$Typ, x)
            r.invscale(x)*(length(r)-1) + 1
        end

        Base.getindex(p::$Typ, i::Int) = p.vals[i]
        Base.length(p::$Typ)  = p.N
        Base.show(io::IO, p::$Typ) = begin
            Typ = $Typ
            print(io, "$Typ($(p.a), $(p.b), $(p.θ), $(p.N))")
        end

        function Base.searchsortedfirst(p::$Typ, x, args...)
            searchsortedfirst(p.vals, x, args...)
        end


    end end

    #-------------------------------------------------------------------------------
    # Exponential Range
    for Typ in [:ExpRange] @eval begin
        struct $Typ{T,TF} <: AbstractRange{T}
            a::T
            b::T
            θ::TF
            N::Int
            vals::Vector{Float64}
            scale::$(scales[Typ]){T,TF}
            invscale::inv($(scales[Typ])){T,TF}
        end

        function $Typ(a,b; θ = 1, N = 100)
            vals    = domscale($Typ)(a,b,θ).(LinRange(0, 1, N))
            scale   = $(scales[Typ])(a, b, θ)
            invscale= inv($(scales[Typ])(a, b, θ))
            T       = reduce(promote_type, map(typeof, (a,b,θ,scale.s)))
            $Typ(convert(T, a), convert(T, b), θ, N, vals, scale, invscale)
        end

        domscale(p::$Typ)      = p.scale
        domscale(::Type{$Typ}) = $(scales[Typ])

        function rangescale(r::$Typ, x)
            r.invscale(x)*(length(r)-1) + 1
        end

        Base.getindex(p::$Typ, i::Int) = p.vals[i]
        Base.length(p::$Typ)  = p.N
        Base.show(io::IO, p::$Typ) = begin
            Typ = $Typ
            print(io, "$Typ($(p.a), $(p.b), $(p.θ), $(p.N))")
        end

        function Base.searchsortedfirst(p::$Typ, x, args...)
            searchsortedfirst(p.vals, x, args...)
        end
    end end

    #-------------------------------------------------------------------------------
    # Finalize Ranges
    function rangescale(r::AbstractRange, x)
        ((x - first(r)))/(last(r) - first(r)) * (length(r) -1) + 1
    end

    rangescale(r::UnitRange, x::Int) = x


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Scaled Interpolation Wraper
    mutable struct ScaledInterpolations{T, N, IT <: Interpolations.AbstractInterpolation,
                                RT <: NTuple{N, AbstractRange},
                                VT <: AbstractArray{T, N}, AT}
        r::RT       # Interpolation grid (ranges)
        v::VT       # Interpolation values
        itp::IT     # Interpolation object
        args::AT    # Arguments to pass to construct interpolant
    end

    function ScaledInterpolations(r, v, args...)
        itp = extrapolate(interpolate(v, args...), Interpolations.Flat())
        return ScaledInterpolations(tuplefy(r), v, itp, args)
    end

    tuplefy(x::Tuple) = x
    tuplefy(x)        = tuple(x)

    """
    ```
    fit!(s::ScaledInterpolations, [v])
    ```
    Re-fit the interpolation after the underlying data has updated.  Works well if
    you use a view into the original array.  If an array of values is passed, copies
    the values to the scaled interpolations value array and then refits the
    interpolation.
    """
    function fit!(s::ScaledInterpolations)
        s.itp = extrapolate(interpolate(s.v, s.args...), Interpolations.Flat())
    end

    function fit!(s::Array{T,N}) where {T <: ScaledInterpolations, N}
        for v in s
            fit!(v)
        end
    end

    dim(::Type{ScaledInterpolations{T, N, IT, RT, VT, AT}}) where {T, N, IT, RT, VT, AT} = N
    Base.eltype(::Type{ScaledInterpolations{T, N, IT, RT, VT, AT}}) where {T, N, IT, RT, VT, AT} = T
    dim(itp::ScaledInterpolations) = dim(typeof(itp))
    Base.eltype(itp::ScaledInterpolations) = eltype(typeof(itp))

    @generated function (s::ScaledInterpolations)(x::Vararg)
        Ns = dim(s)
        Nx = length(x)
        Ns == Nx || begin
            return quote
                Ns = $Ns
                Nx = $Nx
                throw(ArgumentError("Must have $Ns arguments.  You passed $Nx"))
            end
        end

        ex = Expr(:call, :(s.itp))
        for i = 1:Ns
            push!(ex.args, :(rangescale(s.r[$i], x[$i])))
        end
        ex
    end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Scaled Interpolation Shortcuts for Itp_Type
    const CubicSpline    = BSpline(Cubic(Natural(OnGrid())))
    const MonotoneSpline = FritschButlandMonotonicInterpolation()
    const LinearInterp   = BSpline(Linear())
    const LinearInterp2  = BSpline(Linear())

    function ScaledInterpolations(grid, values, type::FritschButlandMonotonicInterpolation)
        return interpolate(collect(grid), values, type)
    end
