const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const grams_per_kilogram = common.grams_per_kilogram;
const pounds_mass_per_kilogram = common.pounds_mass_per_kilogram;

pub const Kg = UnitValueType();
pub const G = UnitValueType();
pub const Lbsm = UnitValueType();

const mass_labels = [_]UnitLabel{
    .{ .abbreviation = "kg", .plural = "kilograms" },
    .{ .abbreviation = "g", .plural = "grams" },
    .{ .abbreviation = "Lbsₘ", .plural = "pounds mass" },
};

pub const Mass = union(enum) {
    kg: Kg,
    g: G,
    lbsm: Lbsm,

    pub fn convertToSiUnit(self: Mass) Kg {
        return switch (self) {
            .kg => |val| val,
            .g => |val| Kg.init(val.value * grams_per_kilogram),
            .lbsm => |val| Kg.init(val.value * pounds_mass_per_kilogram),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &mass_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return mass_labels[0];
    }

    pub fn getValue(self: Mass) f64 {
        return switch (self) {
            .kg => |val| val.value,
            .g => |val| val.value,
            .lbsm => |val| val.value,
        };
    }

    pub fn tryConvert(self: Mass, unit_display: []const u8) ParseUnitError!Mass {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "kg")) return .{ .kg = Kg.init(value_si) };
        if (isMatch(unit_display, "g")) return .{ .g = G.init(value_si / grams_per_kilogram) };
        if (isMatch(unit_display, "Lbsₘ")) return .{ .lbsm = Lbsm.init(value_si / pounds_mass_per_kilogram) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Mass {
        if (isMatch(raw.unit_display, "kg")) return .{ .kg = Kg.init(raw.value) };
        if (isMatch(raw.unit_display, "g")) return .{ .g = G.init(raw.value) };
        if (isMatch(raw.unit_display, "Lbsₘ")) return .{ .lbsm = Lbsm.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Mass) RawUnit {
        return switch (self) {
            .kg => |val| .{ .value = val.value, .unit_display = "kg" },
            .g => |val| .{ .value = val.value, .unit_display = "g" },
            .lbsm => |val| .{ .value = val.value, .unit_display = "Lbsₘ" },
        };
    }
};
