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

using Blake2
using Test

### BEGIN RFC API ###

function selftest_seq(out::Vector{UInt8}, len::Integer, seed::UInt32)
    t::UInt32 = 0
    a::UInt32 = 0xDEAD4BAD * seed           # prime
    b::UInt32 = 1

    if len > 0
        for i = 1:len         # fill the buf
            t = a + b
            a = b
            b = t
            out[i] = (t >> 24) & 0xFF
        end
    end
end

@testset "Blake2b" begin
    # grand hash of hash results
    blake2b_res::Vector{UInt8} = [
        0xC2, 0x3A, 0x78, 0x00, 0xD9, 0x81, 0x23, 0xBD,
        0x10, 0xF5, 0x06, 0xC6, 0x1E, 0x29, 0xDA, 0x56,
        0x03, 0xD7, 0x63, 0xB8, 0xBB, 0xAD, 0x2E, 0x73,
        0x7F, 0x5E, 0x76, 0x5A, 0x7B, 0xCC, 0xD4, 0x75
    ]

    # parameter sets
    b2b_md_len::Vector{UInt64} = [20, 32, 48, 64]
    b2b_in_len::Vector{UInt64} = [0, 3, 128, 129, 255, 1024]

    # 256-bit hash for testing
    md = zeros(UInt8, 64)
    key = zeros(UInt8, 64)
    in = zeros(UInt8, 1024)

    ctx = Blake2.Blake2bInit(32, key, 0)

    for outlen::UInt64 in b2b_md_len
        for inlen::UInt64 in b2b_in_len
            selftest_seq(in, inlen, UInt32(inlen))    # unkeyed, hash
            Blake2.Blake2b!(md, outlen, key, 0, in, inlen)
            Blake2.Blake2bUpdate!(ctx, md, outlen)   # hash the hash

            selftest_seq(key, outlen, UInt32(outlen))  # keyed, hash
            Blake2.Blake2b!(md, outlen, key, outlen, in, inlen)
            Blake2.Blake2bUpdate!(ctx, md, outlen)   # hash the hash
        end
    end

    # compute and compare the hash of hashes
    Blake2.Blake2bFinal!(ctx, md);

    @test md[1:32] == blake2b_res
end

@testset "Blake2s" begin
    # grand hash of hash results
    blake2s_res::Vector{UInt8} = [
        0x6A, 0x41, 0x1F, 0x08, 0xCE, 0x25, 0xAD, 0xCD,
        0xFB, 0x02, 0xAB, 0xA6, 0x41, 0x45, 0x1C, 0xEC,
        0x53, 0xC5, 0x98, 0xB2, 0x4F, 0x4F, 0xC7, 0x87,
        0xFB, 0xDC, 0x88, 0x79, 0x7F, 0x4C, 0x1D, 0xFE
    ]

    # parameter sets
    b2s_md_len::Vector{UInt64} = [16, 20, 28, 32]
    b2s_in_len::Vector{UInt64} = [0, 3, 64, 65, 255, 1024]

    # 256-bit hash for testing
    md = zeros(UInt8, 32)
    key = zeros(UInt8, 32)
    in = zeros(UInt8, 1024)

    ctx = Blake2.Blake2sInit(32, key, 0)

    for outlen::UInt64 in b2s_md_len
        for inlen::UInt64 in b2s_in_len
            selftest_seq(in, inlen, UInt32(inlen))    # unkeyed, hash
            Blake2.Blake2s!(md, outlen, key, 0, in, inlen)
            Blake2.Blake2sUpdate!(ctx, md, outlen)   # hash the hash

            selftest_seq(key, outlen, UInt32(outlen))  # keyed, hash
            Blake2.Blake2s!(md, outlen, key, outlen, in, inlen)
            Blake2.Blake2sUpdate!(ctx, md, outlen)   # hash the hash
        end
    end

    # compute and compare the hash of hashes
    Blake2.Blake2sFinal!(ctx, md);

    @test md[1:32] == blake2s_res
end

### END RFC API ###

@testset "Blake2bABC" begin
    md = zeros(UInt8, 64)
    key = zeros(UInt8, 0)
    in::Vector{UInt8} = [0x61, 0x62, 0x63]

    exp::Vector{UInt8} = [0xBA, 0x80, 0xA5, 0x3F, 0x98, 0x1C, 0x4D, 0x0D, 0x6A, 0x27, 0x97, 0xB6, 0x9F, 0x12, 0xF6, 0xE9,
    0x4C, 0x21, 0x2F, 0x14, 0x68, 0x5A, 0xC4, 0xB7, 0x4B, 0x12, 0xBB, 0x6F, 0xDB, 0xFF, 0xA2, 0xD1,
    0x7D, 0x87, 0xC5, 0x39, 0x2A, 0xAB, 0x79, 0x2D, 0xC2, 0x52, 0xD5, 0xDE, 0x45, 0x33, 0xCC, 0x95,
    0x18, 0xD3, 0x8A, 0xA8, 0xDB, 0xF1, 0x92, 0x5A, 0xB9, 0x23, 0x86, 0xED, 0xD4, 0x00, 0x99, 0x23]

    out = Blake2.Blake2b!(md, 64, key, 0, in, length(in))
    @test out == exp
end