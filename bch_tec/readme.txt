Low Complexity Three-Error-Correcting BCH Decoder

- Decoder utilizes Peterson algorithm computing the error locator polynomial directly (instead of iterative Berlekamp algorithm)
- Instead of Chien search, the decoder utilizes lookup tables (LUT) that stores the roots of the error locator polynomial
- Single and Double error correcting decoders are also supported
- Shortened BCH codes are supported