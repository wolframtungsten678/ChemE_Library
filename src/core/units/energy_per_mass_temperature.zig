const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const absFloat = common.absFloat;
const energy_per_mass_temp_conversion = common.energy_per_mass_temp_conversion;

pub const JPerKgK = UnitValueType();
pub const BtuPerLbsmR = UnitValueType();

const labels = [_]UnitLabel{
    .{ .abbreviation = "J/(kg · K)", .plural = "joules per kilogram kelvin" },
    .{ .abbreviation = "BTU/(Lbsₘ · °R)", .plural = "british thermal units per pounds mass degrees rankine" },
};

pub const EnergyPerMassTemperature = union(enum) {
    j_per_kg_k: JPerKgK,
    btu_per_lbsm_r: BtuPerLbsmR,

    pub fn abs(self: EnergyPerMassTemperature) f64 {
        return absFloat(self.convertToSiUnit().value);
    }

    pub fn sub(self: EnergyPerMassTemperature, other: EnergyPerMassTemperature) EnergyPerMassTemperature {
        const diff = self.convertToSiUnit().sub(other.convertToSiUnit());
        return .{ .j_per_kg_k = diff };
    }

    pub fn convertToSiUnit(self: EnergyPerMassTemperature) JPerKgK {
        return switch (self) {
            .j_per_kg_k => |val| val,
            .btu_per_lbsm_r => |val| JPerKgK.init(val.value * energy_per_mass_temp_conversion),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return labels[0];
    }

    pub fn getValue(self: EnergyPerMassTemperature) f64 {
        return switch (self) {
            .j_per_kg_k => |val| val.value,
            .btu_per_lbsm_r => |val| val.value,
        };
    }

    pub fn tryConvert(self: EnergyPerMassTemperature, unit_display: []const u8) ParseUnitError!EnergyPerMassTemperature {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "J/(kg · K)")) return .{ .j_per_kg_k = JPerKgK.init(value_si) };
        if (isMatch(unit_display, "BTU/(Lbsₘ · °R)")) return .{ .btu_per_lbsm_r = BtuPerLbsmR.init(value_si / energy_per_mass_temp_conversion) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!EnergyPerMassTemperature {
        if (isMatch(raw.unit_display, "J/(kg · K)")) return .{ .j_per_kg_k = JPerKgK.init(raw.value) };
        if (isMatch(raw.unit_display, "BTU/(Lbsₘ · °R)")) return .{ .btu_per_lbsm_r = BtuPerLbsmR.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: EnergyPerMassTemperature) RawUnit {
        return switch (self) {
            .j_per_kg_k => |val| .{ .value = val.value, .unit_display = "J/(kg · K)" },
            .btu_per_lbsm_r => |val| .{ .value = val.value, .unit_display = "BTU/(Lbsₘ · °R)" },
        };
    }
};
