using Blake2
using Test

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