#64-bit Bit Rotations
@inline ror(x::UInt64, k::Int)::UInt64 = (x >>> (0x3f & k)) | (x << (0x3f & -k))
@inline rol(x::UInt64, k::Int)::UInt64 = (x >>> (0x3f & -k)) | (x << (0x3f & k))

#32-bit Bit Rotations
@inline ror(x::UInt32, k::Int)::UInt32 = (x >>> (0x1f & k)) | (x << (0x1f & -k))
@inline rol(x::UInt32, k::Int)::UInt32 = (x >>> (0x1f & -k)) | (x << (0x1f & k))