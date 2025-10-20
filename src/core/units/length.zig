const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const absFloat = common.absFloat;
const meters_per_kilometer = common.meters_per_kilometer;
const feet_per_meter = common.feet_per_meter;
const inches_per_foot = common.inches_per_foot;

pub const M = UnitValueType();
pub const Km = UnitValueType();
pub const Ft = UnitValueType();
pub const Inches = UnitValueType();

const length_labels = [_]UnitLabel{
    .{ .abbreviation = "m", .plural = "meters" },
    .{ .abbreviation = "km", .plural = "kilometers" },
    .{ .abbreviation = "ft", .plural = "feet" },
    .{ .abbreviation = "in", .plural = "inches" },
};

pub const Length = union(enum) {
    m: M,
    km: Km,
    ft: Ft,
    inches: Inches,

    pub fn abs(self: Length) f64 {
        return absFloat(self.convertToSiUnit().value);
    }

    pub fn sub(self: Length, other: Length) Length {
        const diff = self.convertToSiUnit().sub(other.convertToSiUnit());
        return .{ .m = diff };
    }

    pub fn convertToSiUnit(self: Length) M {
        return switch (self) {
            .m => |val| val,
            .km => |val| M.init(val.value * meters_per_kilometer),
            .ft => |val| M.init(val.value / feet_per_meter),
            .inches => |val| M.init(val.value / (feet_per_meter * inches_per_foot)),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &length_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return length_labels[0];
    }

    pub fn getValue(self: Length) f64 {
        return switch (self) {
            .m => |val| val.value,
            .km => |val| val.value,
            .ft => |val| val.value,
            .inches => |val| val.value,
        };
    }

    pub fn tryConvert(self: Length, unit_display: []const u8) ParseUnitError!Length {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "m")) return .{ .m = M.init(value_si) };
        if (isMatch(unit_display, "km")) return .{ .km = Km.init(value_si / meters_per_kilometer) };
        if (isMatch(unit_display, "ft")) return .{ .ft = Ft.init(value_si * feet_per_meter) };
        if (isMatch(unit_display, "in")) return .{ .inches = Inches.init(value_si * feet_per_meter * inches_per_foot) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Length {
        if (isMatch(raw.unit_display, "m")) return .{ .m = M.init(raw.value) };
        if (isMatch(raw.unit_display, "km")) return .{ .km = Km.init(raw.value) };
        if (isMatch(raw.unit_display, "ft")) return .{ .ft = Ft.init(raw.value) };
        if (isMatch(raw.unit_display, "in")) return .{ .inches = Inches.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Length) RawUnit {
        return switch (self) {
            .m => |val| .{ .value = val.value, .unit_display = "m" },
            .km => |val| .{ .value = val.value, .unit_display = "km" },
            .ft => |val| .{ .value = val.value, .unit_display = "ft" },
            .inches => |val| .{ .value = val.value, .unit_display = "in" },
        };
    }
};
