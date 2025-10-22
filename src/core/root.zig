//! By convention, root.zig is the root source file when making a library.
const std = @import("std");

export fn add(a: i32, b: i32) i32 {
    return a + b;
}

test "basic add functionality" {
    try std.testing.expect(add(3, 7) == 10);
}

// Need to do this do that `zig build test` will actually run the tests!
test {
    _ = @import("units.zig");
    _ = @import("thermo/steam/iapws97.zig");
}
