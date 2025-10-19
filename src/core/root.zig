//! Root Zig file for the ChemE dynamic library.
const std = @import("std");

const gsl = @cImport({
    @cInclude("gsl/gsl_interp.h");
});

export fn add(a: i32, b: i32) i32 {
    return a + b;
}

export fn csplineMinSize() usize {
    const raw = gsl.gsl_interp_type_min_size(gsl.gsl_interp_cspline);
    return std.math.cast(usize, raw) orelse @panic("gsl_interp_type_min_size overflow");
}

test "basic add functionality" {
    try std.testing.expect(add(3, 7) == 10);
}

test "basic csplineMinSize functionality" {
    try std.testing.expect(csplineMinSize() == 10);
}
