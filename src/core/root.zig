//! By convention, root.zig is the root source file when making a library.
const std = @import("std");

// Need to do this do that `zig build test` will actually run the tests!
test {
    _ = @import("units.zig");
    _ = @import("thermo/steam/iapws97.zig");
}
