mutable struct Blake2bContext
    b::Vector{UInt8}     # input buffer; len = 128
    h::Vector{UInt64}    # chained state; len = 8
    t::Vector{UInt64}    # total number of bytes; len = 2
    c::UInt64            # pointer for b[]
    outlen::UInt64       # digest size
end

# Mixing function G from the RFC
function Blake2bMix!(v::Vector{UInt64}, a::Integer, b::Integer, c::Integer, d::Integer, x::UInt64, y::UInt64)
    v[a] = v[a] + v[b] + x
    v[d] = ror(xor(v[d], v[a]), 32)
    v[c] = v[c] + v[d]
    v[b] = ror(xor(v[b], v[c]), 24)
    v[a] = v[a] + v[b] + y
    v[d] = ror(xor(v[d], v[a]), 16)
    v[c] = v[c] + v[d]
    v[b] = ror(xor(v[b], v[c]), 63)
end

# Compression function. "last" flag indicates last block.
function Blake2bCompress!(ctx::Blake2bContext, last::Bool)
    sigma::Matrix{UInt8} = [  
     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
    15  11   5   9  10  16  14   7   2  13   1   3  12   8   6   4
    12   9  13   1   6   3  16  14  11  15   4   7   8   2  10   5
     8  10   4   2  14  13  12  15   3   7   6  11   5   1  16   9
    10   1   6   8   3   5  11  16  15   2  12  13   7   9   4  14
     3  13   7  11   1  12   9   4   5  14   8   6  16  15   2  10
    13   6   2  16  15  14   5  11   1   8   7   4  10   3   9  12
    14  12   8  15  13   2   4  10   6   1  16   5   9   7   3  11
     7  16  15  10  12   4   1   9  13   3  14   8   2   5  11   6
    11   3   9   5   8   7   2   6  16  12  10  15   4  13  14   1
     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
    15  11   5   9  10  16  14   7   2  13   1   3  12   8   6   4
    ]

    v::Vector{UInt64} = [ctx.h; [
        0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
        0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
        0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
        0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
    ]]

    v[13] = xor(v[13], ctx.t[1]) # low 64 bits of offset
    v[14] = xor(v[14], ctx.t[2]) # high 64 bits
    if last
        # last block flag set ?
        v[15] = ~v[15]
    end

    m::Vector{UInt64} = reinterpret(UInt64, ctx.b)

    for i = 1:12         # twelve rounds
        Blake2bMix!(v, 1, 5, 9, 13, m[sigma[i, 1]], m[sigma[i, 2]])
        Blake2bMix!(v, 2, 6, 10, 14, m[sigma[i, 3]], m[sigma[i, 4]])
        Blake2bMix!(v, 3, 7, 11, 15, m[sigma[i, 5]], m[sigma[i, 6]])
        Blake2bMix!(v, 4, 8, 12, 16, m[sigma[i, 7]], m[sigma[i, 8]])
        Blake2bMix!(v, 1, 6, 11, 16, m[sigma[i, 9]], m[sigma[i, 10]])
        Blake2bMix!(v, 2, 7, 12, 13, m[sigma[i, 11]], m[sigma[i, 12]])
        Blake2bMix!(v, 3, 8, 9, 14, m[sigma[i, 13]], m[sigma[i, 14]])
        Blake2bMix!(v, 4, 5, 10, 15, m[sigma[i, 15]], m[sigma[i, 16]])
    end

    for i = 1:8
        ctx.h[i] = xor(ctx.h[i], xor(v[i], v[i+8]))
    end
end

# Initialize the hashing context "ctx" with optional key "key".
# 1 <= outlen <= 64 gives the digest size in bytes.
# Secret key (also <= 64 bytes) is optional (keylen = 0).
function Blake2bInit!(ctx::Blake2bContext, outlen::Integer, key::Vector{UInt8}, keylen::Integer)::Blake2bContext # (keylen=0: no key)
    @assert (outlen != 0 && outlen <= 64 && keylen <= 64) "Invalid Parameters"

    ctx.h = [
        0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
        0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
        0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
        0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
    ]
    ctx.h[1] = xor(xor(xor(ctx.h[1], 0x01010000), (keylen << 8)), outlen)
    ctx.t = [0,0]
    ctx.c = 0
    ctx.outlen = outlen
    ctx.b = zeros(UInt8, 128)

    if keylen > 0
        Blake2bUpdate!(ctx, key, keylen)
        ctx.c = 128                  # at the end
    end

    return ctx
end

# Initialize the hashing context "ctx" with optional key "key".
# 1 <= outlen <= 64 gives the digest size in bytes.
# Secret key (also <= 64 bytes) is optional (keylen = 0).
function Blake2bInit(outlen::Integer, key::Vector{UInt8}, keylen::Integer)::Blake2bContext # (keylen=0: no key)
    @assert (outlen != 0 && outlen <= 64 && keylen <= 64) "Invalid Parameters"

    IV::Vector{UInt64} = [
        0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
        0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
        0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
        0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
    ]

    ctx::Blake2bContext = Blake2bContext(zeros(UInt8, 128), IV, zeros(UInt64, 2), 0, outlen)
    ctx.h[1] = xor(xor(xor(ctx.h[1], 0x01010000), (keylen << 8)), outlen)

    if keylen > 0
        Blake2bUpdate!(ctx, key, keylen)
        ctx.c = 128                  # at the end
    end

    return ctx
end

# Add "inlen" bytes from "in" into the hash.
function Blake2bUpdate!(ctx::Blake2bContext, in::Vector{UInt8}, inlen::Integer) # data bytes
    for i = 1:inlen
        if ctx.c == 128                # buffer full ?
            ctx.t[1] += ctx.c          # add counters
            if (ctx.t[1] < ctx.c)       # carry overflow ?
                ctx.t[2] += 1          # high word
            end
            Blake2bCompress!(ctx, false)   # compress (not last)
            ctx.c = 0                  # counter to zero
        end
        ctx.b[ctx.c+1] = in[i]
        ctx.c += 1
    end
end

function Blake2bFinal!(ctx::Blake2bContext, out::Vector{UInt8})
    ctx.t[1] += ctx.c                # mark last block offset
    if (ctx.t[1] < ctx.c)             # carry overflow
        ctx.t[1] += 1                 # high word
    end

    while (ctx.c < 128)                # fill up with zeros
        ctx.b[ctx.c+1] = 0
        ctx.c += 1
    end
    Blake2bCompress!(ctx, true)           # final block flag = 1
    # little endian convert and store
    for i = 0:ctx.outlen-1
        out[i+1] = (ctx.h[(i>>3)+1] >> (8 * (i & 7))) & 0xFF
    end
end

function Blake2b!(out::Vector{UInt8}, outlen::Integer, key::Vector{UInt8}, keylen::Integer, in::Vector{UInt8}, inlen::Integer)
    ctx = Blake2bInit(outlen, key, keylen)
    Blake2bUpdate!(ctx, in, inlen)
    Blake2bFinal!(ctx, out)
    return out
end