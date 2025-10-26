const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const feet_per_meter = common.feet_per_meter;

pub const MPerSec = UnitValueType();
pub const FtPerSec = UnitValueType();

const velocity_labels = [_]UnitLabel{
    .{ .abbreviation = "m/s", .plural = "meters per second" },
    .{ .abbreviation = "ft/s", .plural = "feet per second" },
};

pub const Velocity = union(enum) {
    m_per_sec: MPerSec,
    ft_per_sec: FtPerSec,

    pub fn convertToSiUnit(self: Velocity) MPerSec {
        return switch (self) {
            .m_per_sec => |val| val,
            .ft_per_sec => |val| MPerSec.init(val.value / feet_per_meter),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &velocity_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return velocity_labels[0];
    }

    pub fn getValue(self: Velocity) f64 {
        return switch (self) {
            .m_per_sec => |val| val.value,
            .ft_per_sec => |val| val.value,
        };
    }

    pub fn tryConvert(self: Velocity, unit_display: []const u8) ParseUnitError!Velocity {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "m/s")) return .{ .m_per_sec = MPerSec.init(value_si) };
        if (isMatch(unit_display, "ft/s")) return .{ .ft_per_sec = FtPerSec.init(value_si * feet_per_meter) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Velocity {
        if (isMatch(raw.unit_display, "m/s")) return .{ .m_per_sec = MPerSec.init(raw.value) };
        if (isMatch(raw.unit_display, "ft/s")) return .{ .ft_per_sec = FtPerSec.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Velocity) RawUnit {
        return switch (self) {
            .m_per_sec => |val| .{ .value = val.value, .unit_display = "m/s" },
            .ft_per_sec => |val| .{ .value = val.value, .unit_display = "ft/s" },
        };
    }
};
