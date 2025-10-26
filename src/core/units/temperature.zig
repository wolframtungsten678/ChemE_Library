const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const kelvin_offset_celsius = common.kelvin_offset_celsius;

pub const K = UnitValueType();
pub const C = UnitValueType();
pub const F = UnitValueType();
pub const R = UnitValueType();

const temperature_labels = [_]UnitLabel{
    .{ .abbreviation = "K", .plural = "kelvin" },
    .{ .abbreviation = "°C", .plural = "degrees celsius" },
    .{ .abbreviation = "°F", .plural = "degrees fahrenheit" },
    .{ .abbreviation = "°R", .plural = "degrees rankine" },
};

pub const Temperature = union(enum) {
    k: K,
    c: C,
    f: F,
    r: R,

    pub fn convertToSiUnit(self: Temperature) K {
        return switch (self) {
            .k => |val| val,
            .c => |val| K.init(val.value + kelvin_offset_celsius),
            .f => |val| K.init(((val.value - 32.0) * 5.0 / 9.0) + kelvin_offset_celsius),
            .r => |val| K.init(val.value * 5.0 / 9.0),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &temperature_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return temperature_labels[0];
    }

    pub fn getValue(self: Temperature) f64 {
        return switch (self) {
            .k => |val| val.value,
            .c => |val| val.value,
            .f => |val| val.value,
            .r => |val| val.value,
        };
    }

    pub fn tryConvert(self: Temperature, unit_display: []const u8) ParseUnitError!Temperature {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "K")) return .{ .k = K.init(value_si) };
        if (isMatch(unit_display, "°C")) return .{ .c = C.init(value_si - kelvin_offset_celsius) };
        if (isMatch(unit_display, "°F")) return .{ .f = F.init(((value_si - kelvin_offset_celsius) * 9.0 / 5.0) + 32.0) };
        if (isMatch(unit_display, "°R")) return .{ .r = R.init(value_si * 9.0 / 5.0) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Temperature {
        if (isMatch(raw.unit_display, "K")) return .{ .k = K.init(raw.value) };
        if (isMatch(raw.unit_display, "°C")) return .{ .c = C.init(raw.value) };
        if (isMatch(raw.unit_display, "°F")) return .{ .f = F.init(raw.value) };
        if (isMatch(raw.unit_display, "°R")) return .{ .r = R.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Temperature) RawUnit {
        return switch (self) {
            .k => |val| .{ .value = val.value, .unit_display = "K" },
            .c => |val| .{ .value = val.value, .unit_display = "°C" },
            .f => |val| .{ .value = val.value, .unit_display = "°F" },
            .r => |val| .{ .value = val.value, .unit_display = "°R" },
        };
    }
};
