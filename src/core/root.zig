//! By convention, root.zig is the root source file when making a library.
const std = @import("std");
const iapws97 = @import("thermo/steam/iapws97.zig");
const units = @import("units.zig");

// Need to do this do that `zig build test` will actually run the tests!
test {
    _ = units;
    _ = iapws97;
}

export fn add(a: c_int, b: c_int) c_int {
    return a + b;
}

const Add2 = extern struct {
    a: c_int,
    b: c_int,
    result: c_int,
};

export fn add2(a: c_int, b: c_int) Add2 {
    return .{
        .a = a,
        .b = b,
        .result = a + b,
    };
}
