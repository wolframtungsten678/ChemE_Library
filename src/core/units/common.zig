const utils = @import("unit-utils.zig");

pub const RawUnit = struct {
    value: f64,
    unit_display: []const u8,
};

pub const ParseUnitError = error{
    UnknownUnit,
};

pub const UnitLabel = struct {
    abbreviation: []const u8,
    plural: []const u8,
};

pub const isMatch = utils.isMatch;
pub const absFloat = utils.absFloat;
pub const UnitValueType = utils.UnitValueType;

pub const meters_per_kilometer = 1000.0;
pub const feet_per_meter = 3.28084;
pub const inches_per_foot = 12.0;
pub const meters_per_kilometer_sq = meters_per_kilometer * meters_per_kilometer;
pub const meters_per_kilometer_cu = meters_per_kilometer_sq * meters_per_kilometer;
pub const feet_per_meter_sq = feet_per_meter * feet_per_meter;
pub const feet_per_meter_cu = feet_per_meter_sq * feet_per_meter;
pub const inches_per_meter = feet_per_meter * inches_per_foot;
pub const inches_per_meter_sq = inches_per_meter * inches_per_meter;
pub const inches_per_meter_cu = inches_per_meter_sq * inches_per_meter;
pub const seconds_per_minute = 60.0;
pub const grams_per_kilogram = 1000.0;
pub const pounds_mass_per_kilogram = 2.20462;
pub const kelvin_offset_celsius = 273.15;
pub const pascals_per_kpa = 1000.0;
pub const pascals_per_psi = 6894.76;
pub const energy_per_mass_conversion = 2.2 * 1055.06;
pub const energy_per_mass_temp_conversion = 2.2 * 1055.06 * 5.0 / 9.0;
pub const specific_volume_conversion = 2.20462 / 35.3147;
pub const density_conversion = 35.3147 / 2.20462;
