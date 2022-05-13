# This file contains code originally sourced from code components publised in IETF RFC 7693 (https://datatracker.ietf.org/doc/html/rfc7693).
# The code components are licensed under the Revised BSD License, which is included below:
# "Copyright (c) 2022 IETF Trust and the persons identified as authors of the code. 
# All rights reserved. Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#     - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#     - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
#       in the documentation and/or other materials provided with the distribution.
#     - Neither the name of Internet Society, IETF or IETF Trust, nor the names of specific contributors, may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
#
# Any code contained within this project that is considered original work by the author(s) is licensed and distributed under the BSD 3 Clause License.

using StaticArrays

### BEGIN RFC API ###

struct Blake2sContext
    b::Vector{UInt8}        # input buffer; len = 64
    h::MVector{8,UInt32}    # chained state; len = 8
    t::MVector{2,UInt32}    # total number of bytes; len = 2
    c::MVector{1,UInt64}    # pointer for b[]
    outlen::UInt64          # digest size
end

@inline function B2S_G!(v::Vector{UInt32}, a::Integer, b::Integer, c::Integer, d::Integer, x::UInt32, y::UInt32)
    v[a] = v[a] + v[b] + x
    v[d] = bitrotate(xor(v[d], v[a]), -16)
    v[c] = v[c] + v[d]
    v[b] = bitrotate(xor(v[b], v[c]), -12)
    v[a] = v[a] + v[b] + y
    v[d] = bitrotate(xor(v[d], v[a]), -8)
    v[c] = v[c] + v[d]
    v[b] = bitrotate(xor(v[b], v[c]), -7)
end

# Compression function. "last" flag indicates last block.
function Blake2sCompress!(ctx::Blake2sContext, last::Bool)
    sigma::SMatrix{10, 16, UInt8} = [  
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
    ]

    v::Vector{UInt32} = [ctx.h; [
        0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
        0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
    ]]

    v[13] = xor(v[13], ctx.t[1]) # low 64 bits of offset
    v[14] = xor(v[14], ctx.t[2]) # high 64 bits
    if last
        # last block flag set ?
        v[15] = ~v[15]
    end

    m::Vector{UInt32} = reinterpret(UInt32, ctx.b)

    for i = 1:10         # twelve rounds
        B2S_G!(v, 1, 5, 9, 13, m[sigma[i, 1]], m[sigma[i, 2]])
        B2S_G!(v, 2, 6, 10, 14, m[sigma[i, 3]], m[sigma[i, 4]])
        B2S_G!(v, 3, 7, 11, 15, m[sigma[i, 5]], m[sigma[i, 6]])
        B2S_G!(v, 4, 8, 12, 16, m[sigma[i, 7]], m[sigma[i, 8]])
        B2S_G!(v, 1, 6, 11, 16, m[sigma[i, 9]], m[sigma[i, 10]])
        B2S_G!(v, 2, 7, 12, 13, m[sigma[i, 11]], m[sigma[i, 12]])
        B2S_G!(v, 3, 8, 9, 14, m[sigma[i, 13]], m[sigma[i, 14]])
        B2S_G!(v, 4, 5, 10, 15, m[sigma[i, 15]], m[sigma[i, 16]])
    end

    for i = 1:8
        ctx.h[i] = xor(ctx.h[i], xor(v[i], v[i+8]))
    end
end

# Initialize the hashing context "ctx" with optional key "key".
# 1 <= outlen <= 64 gives the digest size in bytes.
# Secret key (also <= 64 bytes) is optional (keylen = 0).
function Blake2sInit!(ctx::Blake2sContext, outlen::Integer, key::Vector{UInt8}, keylen::Integer)::Blake2sContext # (keylen=0: no key)
    @assert (outlen != 0 && outlen <= 32 && keylen <= 32) "Invalid Parameters"

    ctx.h = [0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A, 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19]
    ctx.h[1] = xor(xor(xor(ctx.h[1], 0x01010000), (keylen << 8)), outlen)
    ctx.t = [0,0]
    ctx.c = [0]
    ctx.outlen = outlen
    ctx.b = zeros(UInt8, 64)

    if keylen > 0
        Blake2sUpdate!(ctx, key, keylen)
        ctx.c[1] = 64                  # at the end
    end

    return ctx
end

# Add "inlen" bytes from "in" into the hash.
function Blake2sUpdate!(ctx::Blake2sContext, in::Vector{UInt8}, inlen::Integer) # data bytes
    for i = 1:inlen
        if ctx.c[1] == 64                # buffer full ?
            ctx.t[1] += ctx.c[1]          # add counters
            if (ctx.t[1] < ctx.c[1])       # carry overflow ?
                ctx.t[2] += 1          # high word
            end
            Blake2sCompress!(ctx, false)   # compress (not last)
            ctx.c[1] = 0                  # counter to zero
        end
        ctx.b[ctx.c[1]+1] = in[i]
        ctx.c[1] += 1
    end
end

function Blake2sFinal!(ctx::Blake2sContext, out::Vector{UInt8})
    ctx.t[1] += ctx.c[1]                # mark last block offset
    if (ctx.t[1] < ctx.c[1])             # carry overflow
        ctx.t[1] += 1                 # high word
    end

    while (ctx.c[1] < 64)                # fill up with zeros
        ctx.b[ctx.c[1]+1] = 0
        ctx.c[1] += 1
    end
    Blake2sCompress!(ctx, true)           # final block flag = 1
    # little endian convert and store
    for i = 0:ctx.outlen-1
        out[i+1] = (ctx.h[(i>>2)+1] >> (8 * (i & 3))) & 0xFF
    end
end

function Blake2s!(out::Vector{UInt8}, outlen::Integer, key::Vector{UInt8}, keylen::Integer, in::Vector{UInt8}, inlen::Integer)
    ctx = Blake2sInit(outlen, key, keylen)
    Blake2sUpdate!(ctx, in, inlen)
    Blake2sFinal!(ctx, out)
    return out
end

### END RFC API ###

# Variant on the above init function which returns a context rather modifying a passed context.
function Blake2sInit(outlen::Integer, key::Vector{UInt8}, keylen::Integer)::Blake2sContext # (keylen=0: no key)
    @assert (outlen != 0 && outlen <= 32 && keylen <= 32) "Invalid Parameters"

    IV::MVector{8, UInt32} = [
        0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
        0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
    ]

    ctx::Blake2sContext = Blake2sContext(zeros(UInt8, 64), IV, zeros(UInt32, 2), zeros(UInt32, 1), outlen)
    ctx.h[1] = xor(xor(xor(ctx.h[1], 0x01010000), (keylen << 8)), outlen)

    if keylen > 0
        Blake2sUpdate!(ctx, key, keylen)
        ctx.c[1] = 64                  # at the end
    end

    return ctx
end