#pragma once

namespace popvcf
{
//! Buffer size when encoding
long constexpr ENC_BUFFER_SIZE{65536};

//! Writes the encoding of a stdin VCF to stdout
void encode();

} // namespace popvcf
