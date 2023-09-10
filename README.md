This is IEEE APSqRt, a square root implementation for use with [rustc\_apfloat](https://crates.io/crates/rustc_apfloat)'s IEEE 754 floats.

(APFloat is short for Arbitrary Precision Floating Point. IEEE APSqRt is short for IEEE 754 APFloat Square Root. Unlike APFloat, it is only designed to work with IEEE 754 single-, double-, and quad-precision floats; thus, it does not provide arbitrary precision itself.)

# Why

rustc\_apfloat is a software floating point library. It's a Rust conversion of "APFloat" from the LLVM project. These libraries emphasize clarity and correctness over speed, but are still pretty fast. They're used, respectively, in the `rustc` and `clang` compilers to do compile-time calculation of floating point operations. These compilers don't use the host's floating point because, when cross compiling, host and target behavior may be difficult or impossible to reconcile. Emulators often have the same problem. I was writing a deterministic RV32-G emulator, and as part of that, I needed to emulate the RISC-V standard's floating point specifications specifically. rustc\_apfloat to the rescue!

Because rustc\_apfloat implements almost every operation needed for RISC-V floating point emulation, it made my life a lot easier—but that last operation, square root, is important, so I had to implement it myself. It would have been faster and easier to use the host's `sqrt` functions for these, and to accept the inconsistencies in emulation, but my emulator is meant to fit in a multiplayer game's logic loop (?!!?!) so it *has* to be deterministic. A bad algorithm that is wrong in exactly the same ways on any platform is, therefore, better than a good algorithm that has even the slightest disagreement between platforms.

I originally implemented a very naive guess and check algorithm (see [Bad APSqRt][2]), but once I'd done that, I had everything I needed to implement and understand Newton-Raphson instead. This crate is the result.

[2]: https://github.com/SolraBizna/bad-apsqrt

# Algorithm

IEEE APSqRt implements the Newton-Raphson method, a.k.a. the Babylonian method. It forms an initial guess of `2^(exponent/2)`, and repeatedly refines it using the formula: `new_guess = (old_guess + square / old_guess) / 2`. After just 5-7 iterations, this gets as close as it's going to get to an answer.

# Accuracy

If there is an exact solution, either form of IEEE APSqRt always finds it. For `sqrt_fast`, roughly 70% of inexact solutions will be correct, roughly 30% of solutions will be off by one ulp, and a tiny percentage of solutions will be off by two ulps. For `sqrt_slow`, which performs the operations at higher precision, it will usually take only one more iteration and the answer will be 100% correct. (`sqrt_slow` is, unfortunately, due to limitations of rustc\_apfloat's external interface, only available for 32- and 64-bit floats.)

All NaNs produced by IEEE APSqRt are "canon NaNs" according to the RISC-V standard.

# How to use

Call `ieee_apsqrt::sqrt_fast` with a `u32`/`u64`/`u128`, or `ieee_apsqrt::sqrt_accurate` with a `u32`/`u64`. The former function is a little less accurate, but about twice as fast. The latter function is as accurate as possible, but about half as fast. Both functions also require a rounding mode.

Both functions return a tuple of `(rustc_apfloat::StatusAnd<uXX>, u32)`, where the first value is the result of the calculation, and the second value is how many iterations of Newton-Raphson were required.

If you're counting clock cycles for emulation purposes, consider the square root to have a base cost of a single multiply, and an iteration cost of one division, one addition, and one multiply. Consider `sqrt_fast` to be performing operations at the requested precision, and `sqrt_accurate` at twice that precision (e.g. single -> double, double -> quad).

# Legalese

IEEE APSqRt is copyright 2023 Solra Bizna.

IEEE APSqRt is licensed under the Apache 2 with LLVM exception license. This is the same license as the rustc\_apfloat crate. This is the simplest way to guarantee that IEEE APSqRt can be used anywhere anywhere rustc\_apfloat can.
