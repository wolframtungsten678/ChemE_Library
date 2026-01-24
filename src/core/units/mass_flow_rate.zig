const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const seconds_per_minute = common.seconds_per_minute;
const pounds_mass_per_kilogram = common.pounds_mass_per_kilogram;

pub const KgPerSec = UnitValueType();
pub const KgPerMin = UnitValueType();
pub const LbsmPerSec = UnitValueType();
pub const LbsmPerMin = UnitValueType();

const mass_flow_rate_labels = [_]UnitLabel{
    .{ .abbreviation = "kg/sec", .plural = "kilograms per second" },
    .{ .abbreviation = "kg/min", .plural = "kilograms per minute" },
    .{ .abbreviation = "lbsm/sec", .plural = "pounds mass per second" },
    .{ .abbreviation = "lbsm/min", .plural = "pounds mass per minute" },
};

pub const MassFlowRate = union(enum) {
    kg_per_sec: KgPerSec,
    kg_per_min: KgPerMin,
    lbsm_per_sec: LbsmPerSec,
    lbsm_per_min: LbsmPerMin,

    pub fn convertToSiUnit(self: MassFlowRate) KgPerSec {
        return switch (self) {
            .kg_per_sec => |val| val,
            .kg_per_min => |val| KgPerSec.init(val.value / seconds_per_minute),
            .lbsm_per_sec => |val| KgPerSec.init(val.value / pounds_mass_per_kilogram),
            .lbsm_per_min => |val| KgPerSec.init(val.value / (pounds_mass_per_kilogram * seconds_per_minute)),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &mass_flow_rate_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return mass_flow_rate_labels[0];
    }

    pub fn getValue(self: MassFlowRate) f64 {
        return switch (self) {
            .kg_per_sec => |val| val.value,
            .kg_per_min => |val| val.value,
            .lbsm_per_sec => |val| val.value,
            .lbsm_per_min => |val| val.value,
        };
    }

    pub fn tryConvert(self: MassFlowRate, unit_display: []const u8) ParseUnitError!MassFlowRate {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "kg/sec")) return .{ .kg_per_sec = KgPerSec.init(value_si) };
        if (isMatch(unit_display, "kg/min")) return .{ .kg_per_min = KgPerMin.init(value_si * seconds_per_minute) };
        if (isMatch(unit_display, "lbsm/sec")) return .{ .lbsm_per_sec = LbsmPerSec.init(value_si * pounds_mass_per_kilogram) };
        if (isMatch(unit_display, "lbsm/min")) return .{ .lbsm_per_min = LbsmPerMin.init(value_si * pounds_mass_per_kilogram * seconds_per_minute) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!MassFlowRate {
        if (isMatch(raw.unit_display, "kg/sec")) return .{ .kg_per_sec = KgPerSec.init(raw.value) };
        if (isMatch(raw.unit_display, "kg/min")) return .{ .kg_per_min = KgPerMin.init(raw.value) };
        if (isMatch(raw.unit_display, "lbsm/sec")) return .{ .lbsm_per_sec = LbsmPerSec.init(raw.value) };
        if (isMatch(raw.unit_display, "lbsm/min")) return .{ .lbsm_per_min = LbsmPerMin.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: MassFlowRate) RawUnit {
        return switch (self) {
            .kg_per_sec => |val| .{ .value = val.value, .unit_display = "kg/sec" },
            .kg_per_min => |val| .{ .value = val.value, .unit_display = "kg/min" },
            .lbsm_per_sec => |val| .{ .value = val.value, .unit_display = "lbsm/sec" },
            .lbsm_per_min => |val| .{ .value = val.value, .unit_display = "lbsm/min" },
        };
    }
};
