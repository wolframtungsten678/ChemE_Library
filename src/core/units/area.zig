const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const meters_per_kilometer_sq = common.meters_per_kilometer_sq;
const feet_per_meter_sq = common.feet_per_meter_sq;
const inches_per_meter_sq = common.inches_per_meter_sq;

pub const M2 = UnitValueType();
pub const Km2 = UnitValueType();
pub const Ft2 = UnitValueType();
pub const Inches2 = UnitValueType();

const area_labels = [_]UnitLabel{
    .{ .abbreviation = "m²", .plural = "meters squared" },
    .{ .abbreviation = "km²", .plural = "kilometers squared" },
    .{ .abbreviation = "ft²", .plural = "feet squared" },
    .{ .abbreviation = "in²", .plural = "inches squared" },
};

pub const Area = union(enum) {
    m2: M2,
    km2: Km2,
    ft2: Ft2,
    inches2: Inches2,

    pub fn convertToSiUnit(self: Area) M2 {
        return switch (self) {
            .m2 => |val| val,
            .km2 => |val| M2.init(val.value * meters_per_kilometer_sq),
            .ft2 => |val| M2.init(val.value / feet_per_meter_sq),
            .inches2 => |val| M2.init(val.value / inches_per_meter_sq),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &area_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return area_labels[0];
    }

    pub fn getValue(self: Area) f64 {
        return switch (self) {
            .m2 => |val| val.value,
            .km2 => |val| val.value,
            .ft2 => |val| val.value,
            .inches2 => |val| val.value,
        };
    }

    pub fn tryConvert(self: Area, unit_display: []const u8) ParseUnitError!Area {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "m²")) return .{ .m2 = M2.init(value_si) };
        if (isMatch(unit_display, "km²")) return .{ .km2 = Km2.init(value_si / meters_per_kilometer_sq) };
        if (isMatch(unit_display, "ft²")) return .{ .ft2 = Ft2.init(value_si * feet_per_meter_sq) };
        if (isMatch(unit_display, "in²")) return .{ .inches2 = Inches2.init(value_si * inches_per_meter_sq) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Area {
        if (isMatch(raw.unit_display, "m²")) return .{ .m2 = M2.init(raw.value) };
        if (isMatch(raw.unit_display, "km²")) return .{ .km2 = Km2.init(raw.value) };
        if (isMatch(raw.unit_display, "ft²")) return .{ .ft2 = Ft2.init(raw.value) };
        if (isMatch(raw.unit_display, "in²")) return .{ .inches2 = Inches2.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Area) RawUnit {
        return switch (self) {
            .m2 => |val| .{ .value = val.value, .unit_display = "m²" },
            .km2 => |val| .{ .value = val.value, .unit_display = "km²" },
            .ft2 => |val| .{ .value = val.value, .unit_display = "ft²" },
            .inches2 => |val| .{ .value = val.value, .unit_display = "in²" },
        };
    }
};
