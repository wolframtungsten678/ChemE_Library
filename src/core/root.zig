//! By convention, root.zig is the root source file when making a library.
const std = @import("std");
const iapws97 = @import("thermo/steam/iapws97.zig");
const units = @import("units.zig");

// Need to do this do that `zig build test` will actually run the tests!
test {
    _ = units;
    _ = iapws97;
}
