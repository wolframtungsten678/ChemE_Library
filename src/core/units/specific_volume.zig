const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const absFloat = common.absFloat;
const specific_volume_conversion = common.specific_volume_conversion;

pub const M3PerKg = UnitValueType();
pub const Ft3PerLbsm = UnitValueType();

const labels = [_]UnitLabel{
    .{ .abbreviation = "m³/kg", .plural = "cubic meters per kilogram" },
    .{ .abbreviation = "ft³/Lbsₘ", .plural = "cubic feet per pounds mass" },
};

pub const SpecificVolume = union(enum) {
    m3_per_kg: M3PerKg,
    ft3_per_lbsm: Ft3PerLbsm,

    pub fn abs(self: SpecificVolume) f64 {
        return absFloat(self.convertToSiUnit().value);
    }

    pub fn sub(self: SpecificVolume, other: SpecificVolume) SpecificVolume {
        const diff = self.convertToSiUnit().sub(other.convertToSiUnit());
        return .{ .m3_per_kg = diff };
    }

    pub fn convertToSiUnit(self: SpecificVolume) M3PerKg {
        return switch (self) {
            .m3_per_kg => |val| val,
            .ft3_per_lbsm => |val| M3PerKg.init(val.value * specific_volume_conversion),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return labels[0];
    }

    pub fn getValue(self: SpecificVolume) f64 {
        return switch (self) {
            .m3_per_kg => |val| val.value,
            .ft3_per_lbsm => |val| val.value,
        };
    }

    pub fn tryConvert(self: SpecificVolume, unit_display: []const u8) ParseUnitError!SpecificVolume {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "m³/kg")) return .{ .m3_per_kg = M3PerKg.init(value_si) };
        if (isMatch(unit_display, "ft³/Lbsₘ")) return .{ .ft3_per_lbsm = Ft3PerLbsm.init(value_si / specific_volume_conversion) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!SpecificVolume {
        if (isMatch(raw.unit_display, "m³/kg")) return .{ .m3_per_kg = M3PerKg.init(raw.value) };
        if (isMatch(raw.unit_display, "ft³/Lbsₘ")) return .{ .ft3_per_lbsm = Ft3PerLbsm.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: SpecificVolume) RawUnit {
        return switch (self) {
            .m3_per_kg => |val| .{ .value = val.value, .unit_display = "m³/kg" },
            .ft3_per_lbsm => |val| .{ .value = val.value, .unit_display = "ft³/Lbsₘ" },
        };
    }
};
