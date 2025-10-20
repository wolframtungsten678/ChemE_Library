const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const absFloat = common.absFloat;
const density_conversion = common.density_conversion;

pub const KgPerM3 = UnitValueType();
pub const LbsmPerFt3 = UnitValueType();

const labels = [_]UnitLabel{
    .{ .abbreviation = "kg/m³", .plural = "cubic meters per kilogram" },
    .{ .abbreviation = "Lbsₘ/ft³", .plural = "pounds mass per cubic feet" },
};

pub const Density = union(enum) {
    kg_per_m3: KgPerM3,
    lbsm_per_ft3: LbsmPerFt3,

    pub fn abs(self: Density) f64 {
        return absFloat(self.convertToSiUnit().value);
    }

    pub fn sub(self: Density, other: Density) Density {
        const diff = self.convertToSiUnit().sub(other.convertToSiUnit());
        return .{ .kg_per_m3 = diff };
    }

    pub fn convertToSiUnit(self: Density) KgPerM3 {
        return switch (self) {
            .kg_per_m3 => |val| val,
            .lbsm_per_ft3 => |val| KgPerM3.init(val.value * density_conversion),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return labels[0];
    }

    pub fn getValue(self: Density) f64 {
        return switch (self) {
            .kg_per_m3 => |val| val.value,
            .lbsm_per_ft3 => |val| val.value,
        };
    }

    pub fn tryConvert(self: Density, unit_display: []const u8) ParseUnitError!Density {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "kg/m³")) return .{ .kg_per_m3 = KgPerM3.init(value_si) };
        if (isMatch(unit_display, "Lbsₘ/ft³")) return .{ .lbsm_per_ft3 = LbsmPerFt3.init(value_si / density_conversion) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Density {
        if (isMatch(raw.unit_display, "kg/m³")) return .{ .kg_per_m3 = KgPerM3.init(raw.value) };
        if (isMatch(raw.unit_display, "Lbsₘ/ft³")) return .{ .lbsm_per_ft3 = LbsmPerFt3.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Density) RawUnit {
        return switch (self) {
            .kg_per_m3 => |val| .{ .value = val.value, .unit_display = "kg/m³" },
            .lbsm_per_ft3 => |val| .{ .value = val.value, .unit_display = "Lbsₘ/ft³" },
        };
    }
};
