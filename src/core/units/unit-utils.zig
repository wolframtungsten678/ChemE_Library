const std = @import("std");

pub fn isMatch(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

pub fn UnitValueType() type {
    return struct {
        value: f64,

        pub fn init(value: f64) @This() {
            return .{ .value = value };
        }
    };
}
