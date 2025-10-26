const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const meters_per_kilometer_cu = common.meters_per_kilometer_cu;
const feet_per_meter_cu = common.feet_per_meter_cu;
const inches_per_meter_cu = common.inches_per_meter_cu;

pub const M3 = UnitValueType();
pub const Km3 = UnitValueType();
pub const Ft3 = UnitValueType();
pub const Inches3 = UnitValueType();

const volume_labels = [_]UnitLabel{
    .{ .abbreviation = "m³", .plural = "meters cubed" },
    .{ .abbreviation = "km³", .plural = "kilometers cubed" },
    .{ .abbreviation = "ft³", .plural = "feet cubed" },
    .{ .abbreviation = "in³", .plural = "inches cubed" },
};

pub const Volume = union(enum) {
    m3: M3,
    km3: Km3,
    ft3: Ft3,
    inches3: Inches3,

    pub fn convertToSiUnit(self: Volume) M3 {
        return switch (self) {
            .m3 => |val| val,
            .km3 => |val| M3.init(val.value * meters_per_kilometer_cu),
            .ft3 => |val| M3.init(val.value / feet_per_meter_cu),
            .inches3 => |val| M3.init(val.value / inches_per_meter_cu),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &volume_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return volume_labels[0];
    }

    pub fn getValue(self: Volume) f64 {
        return switch (self) {
            .m3 => |val| val.value,
            .km3 => |val| val.value,
            .ft3 => |val| val.value,
            .inches3 => |val| val.value,
        };
    }

    pub fn tryConvert(self: Volume, unit_display: []const u8) ParseUnitError!Volume {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "m³")) return .{ .m3 = M3.init(value_si) };
        if (isMatch(unit_display, "km³")) return .{ .km3 = Km3.init(value_si / meters_per_kilometer_cu) };
        if (isMatch(unit_display, "ft³")) return .{ .ft3 = Ft3.init(value_si * feet_per_meter_cu) };
        if (isMatch(unit_display, "in³")) return .{ .inches3 = Inches3.init(value_si * inches_per_meter_cu) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Volume {
        if (isMatch(raw.unit_display, "m³")) return .{ .m3 = M3.init(raw.value) };
        if (isMatch(raw.unit_display, "km³")) return .{ .km3 = Km3.init(raw.value) };
        if (isMatch(raw.unit_display, "ft³")) return .{ .ft3 = Ft3.init(raw.value) };
        if (isMatch(raw.unit_display, "in³")) return .{ .inches3 = Inches3.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Volume) RawUnit {
        return switch (self) {
            .m3 => |val| .{ .value = val.value, .unit_display = "m³" },
            .km3 => |val| .{ .value = val.value, .unit_display = "km³" },
            .ft3 => |val| .{ .value = val.value, .unit_display = "ft³" },
            .inches3 => |val| .{ .value = val.value, .unit_display = "in³" },
        };
    }
};
