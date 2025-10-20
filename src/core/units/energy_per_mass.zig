const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const absFloat = common.absFloat;
const energy_per_mass_conversion = common.energy_per_mass_conversion;

pub const JPerKg = UnitValueType();
pub const BtuPerLbsm = UnitValueType();

const energy_per_mass_labels = [_]UnitLabel{
    .{ .abbreviation = "J/kg", .plural = "joules pr kilogram" },
    .{ .abbreviation = "BTU/Lbsₘ", .plural = "british thermal units per pounds mass" },
};

pub const EnergyPerMass = union(enum) {
    j_per_kg: JPerKg,
    btu_per_lbsm: BtuPerLbsm,

    pub fn abs(self: EnergyPerMass) f64 {
        return absFloat(self.convertToSiUnit().value);
    }

    pub fn sub(self: EnergyPerMass, other: EnergyPerMass) EnergyPerMass {
        const diff = self.convertToSiUnit().sub(other.convertToSiUnit());
        return .{ .j_per_kg = diff };
    }

    pub fn convertToSiUnit(self: EnergyPerMass) JPerKg {
        return switch (self) {
            .j_per_kg => |val| val,
            .btu_per_lbsm => |val| JPerKg.init(val.value * energy_per_mass_conversion),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &energy_per_mass_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return energy_per_mass_labels[0];
    }

    pub fn getValue(self: EnergyPerMass) f64 {
        return switch (self) {
            .j_per_kg => |val| val.value,
            .btu_per_lbsm => |val| val.value,
        };
    }

    pub fn tryConvert(self: EnergyPerMass, unit_display: []const u8) ParseUnitError!EnergyPerMass {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "J/kg")) return .{ .j_per_kg = JPerKg.init(value_si) };
        if (isMatch(unit_display, "BTU/Lbsₘ")) return .{ .btu_per_lbsm = BtuPerLbsm.init(value_si / energy_per_mass_conversion) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!EnergyPerMass {
        if (isMatch(raw.unit_display, "J/kg")) return .{ .j_per_kg = JPerKg.init(raw.value) };
        if (isMatch(raw.unit_display, "BTU/Lbsₘ")) return .{ .btu_per_lbsm = BtuPerLbsm.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: EnergyPerMass) RawUnit {
        return switch (self) {
            .j_per_kg => |val| .{ .value = val.value, .unit_display = "J/kg" },
            .btu_per_lbsm => |val| .{ .value = val.value, .unit_display = "BTU/Lbsₘ" },
        };
    }
};
