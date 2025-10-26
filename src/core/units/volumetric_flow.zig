const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const feet_per_meter_cu = common.feet_per_meter_cu;
const seconds_per_minute = common.seconds_per_minute;

pub const M3PerSec = UnitValueType();
pub const M3PerMin = UnitValueType();
pub const Ft3PerSec = UnitValueType();
pub const Ft3PerMin = UnitValueType();

const volumetric_flow_labels = [_]UnitLabel{
    .{ .abbreviation = "m³/sec", .plural = "meters cubed per second" },
    .{ .abbreviation = "m³/min", .plural = "meters cubed per minute" },
    .{ .abbreviation = "ft³/sec", .plural = "feet cubed per second" },
    .{ .abbreviation = "ft³/min", .plural = "feet cubed per minute" },
};

pub const VolumetricFlowRate = union(enum) {
    m3_per_sec: M3PerSec,
    m3_per_min: M3PerMin,
    ft3_per_sec: Ft3PerSec,
    ft3_per_min: Ft3PerMin,

    pub fn convertToSiUnit(self: VolumetricFlowRate) M3PerSec {
        return switch (self) {
            .m3_per_sec => |val| val,
            .m3_per_min => |val| M3PerSec.init(val.value * seconds_per_minute),
            .ft3_per_sec => |val| M3PerSec.init(val.value / feet_per_meter_cu),
            .ft3_per_min => |val| M3PerSec.init(val.value * (seconds_per_minute / feet_per_meter_cu)),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &volumetric_flow_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return volumetric_flow_labels[0];
    }

    pub fn getValue(self: VolumetricFlowRate) f64 {
        return switch (self) {
            .m3_per_sec => |val| val.value,
            .m3_per_min => |val| val.value,
            .ft3_per_sec => |val| val.value,
            .ft3_per_min => |val| val.value,
        };
    }

    pub fn tryConvert(self: VolumetricFlowRate, unit_display: []const u8) ParseUnitError!VolumetricFlowRate {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "m³/sec")) return .{ .m3_per_sec = M3PerSec.init(value_si) };
        if (isMatch(unit_display, "m³/min")) return .{ .m3_per_min = M3PerMin.init(value_si / seconds_per_minute) };
        if (isMatch(unit_display, "ft³/sec")) return .{ .ft3_per_sec = Ft3PerSec.init(value_si * feet_per_meter_cu) };
        if (isMatch(unit_display, "ft³/min")) return .{ .ft3_per_min = Ft3PerMin.init(value_si * (feet_per_meter_cu / seconds_per_minute)) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!VolumetricFlowRate {
        if (isMatch(raw.unit_display, "m³/sec")) return .{ .m3_per_sec = M3PerSec.init(raw.value) };
        if (isMatch(raw.unit_display, "m³/min")) return .{ .m3_per_min = M3PerMin.init(raw.value) };
        if (isMatch(raw.unit_display, "ft³/sec")) return .{ .ft3_per_sec = Ft3PerSec.init(raw.value) };
        if (isMatch(raw.unit_display, "ft³/min")) return .{ .ft3_per_min = Ft3PerMin.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: VolumetricFlowRate) RawUnit {
        return switch (self) {
            .m3_per_sec => |val| .{ .value = val.value, .unit_display = "m³/sec" },
            .m3_per_min => |val| .{ .value = val.value, .unit_display = "m³/min" },
            .ft3_per_sec => |val| .{ .value = val.value, .unit_display = "ft³/sec" },
            .ft3_per_min => |val| .{ .value = val.value, .unit_display = "ft³/min" },
        };
    }
};
