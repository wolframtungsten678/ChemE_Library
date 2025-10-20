const std = @import("std");

pub fn isMatch(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

pub fn absFloat(value: f64) f64 {
    return std.math.fabs(value);
}

pub fn UnitValueType() type {
    return struct {
        value: f64,

        pub fn init(value: f64) @This() {
            return .{ .value = value };
        }

        pub fn abs(self: @This()) f64 {
            return absFloat(self.value);
        }

        pub fn sub(self: @This(), other: @This()) @This() {
            return .{ .value = self.value - other.value };
        }
    };
}
