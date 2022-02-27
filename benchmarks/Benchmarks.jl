module Blake2Benchmarks

using BenchmarkTools
using Blake2

suite = BenchmarkGroup()
suite["Blake2b"] = BenchmarkGroup(["Blake2b512Keyless","Blake2b512Keyed"])

### Blake2b512Keyless
md = zeros(UInt8, 64)
key = zeros(UInt8, 0)
in = rand(UInt8, 64)

suite["Blake2b"]["Blake2b512Keyless"] = @benchmarkable Blake2.Blake2b!($md, length($md), $key, length($key), $in, length($in))

# Run All Benchmarks

function BenchmarkAll()
    tune!(suite)
end

end