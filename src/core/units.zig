const std = @import("std");

const common = @import("units/common.zig");
const length_mod = @import("units/length.zig");
const area_mod = @import("units/area.zig");
const volume_mod = @import("units/volume.zig");
const volumetric_flow_mod = @import("units/volumetric_flow.zig");
const mass_mod = @import("units/mass.zig");
const temperature_mod = @import("units/temperature.zig");
const pressure_mod = @import("units/pressure.zig");
const energy_per_mass_mod = @import("units/energy_per_mass.zig");
const energy_per_mass_temperature_mod = @import("units/energy_per_mass_temperature.zig");
const velocity_mod = @import("units/velocity.zig");
const specific_volume_mod = @import("units/specific_volume.zig");
const density_mod = @import("units/density.zig");

pub const RawUnit = common.RawUnit;
pub const ParseUnitError = common.ParseUnitError;
pub const UnitLabel = common.UnitLabel;
pub const meters_per_kilometer = common.meters_per_kilometer;
pub const feet_per_meter = common.feet_per_meter;
pub const inches_per_foot = common.inches_per_foot;
pub const meters_per_kilometer_sq = common.meters_per_kilometer_sq;
pub const meters_per_kilometer_cu = common.meters_per_kilometer_cu;
pub const feet_per_meter_sq = common.feet_per_meter_sq;
pub const feet_per_meter_cu = common.feet_per_meter_cu;
pub const inches_per_meter = common.inches_per_meter;
pub const inches_per_meter_sq = common.inches_per_meter_sq;
pub const inches_per_meter_cu = common.inches_per_meter_cu;
pub const seconds_per_minute = common.seconds_per_minute;
pub const grams_per_kilogram = common.grams_per_kilogram;
pub const pounds_mass_per_kilogram = common.pounds_mass_per_kilogram;
pub const kelvin_offset_celsius = common.kelvin_offset_celsius;
pub const pascals_per_kpa = common.pascals_per_kpa;
pub const pascals_per_psi = common.pascals_per_psi;
pub const energy_per_mass_conversion = common.energy_per_mass_conversion;
pub const energy_per_mass_temp_conversion = common.energy_per_mass_temp_conversion;
pub const specific_volume_conversion = common.specific_volume_conversion;
pub const density_conversion = common.density_conversion;

pub const M = length_mod.M;
pub const Km = length_mod.Km;
pub const Ft = length_mod.Ft;
pub const Inches = length_mod.Inches;
pub const Length = length_mod.Length;

pub const M2 = area_mod.M2;
pub const Km2 = area_mod.Km2;
pub const Ft2 = area_mod.Ft2;
pub const Inches2 = area_mod.Inches2;
pub const Area = area_mod.Area;

pub const M3 = volume_mod.M3;
pub const Km3 = volume_mod.Km3;
pub const Ft3 = volume_mod.Ft3;
pub const Inches3 = volume_mod.Inches3;
pub const Volume = volume_mod.Volume;

pub const M3PerSec = volumetric_flow_mod.M3PerSec;
pub const M3PerMin = volumetric_flow_mod.M3PerMin;
pub const Ft3PerSec = volumetric_flow_mod.Ft3PerSec;
pub const Ft3PerMin = volumetric_flow_mod.Ft3PerMin;
pub const VolumetricFlowRate = volumetric_flow_mod.VolumetricFlowRate;

pub const Kg = mass_mod.Kg;
pub const G = mass_mod.G;
pub const Lbsm = mass_mod.Lbsm;
pub const Mass = mass_mod.Mass;

pub const K = temperature_mod.K;
pub const C = temperature_mod.C;
pub const F = temperature_mod.F;
pub const R = temperature_mod.R;
pub const Temperature = temperature_mod.Temperature;

pub const Pa = pressure_mod.Pa;
pub const KPa = pressure_mod.KPa;
pub const Lbf = pressure_mod.Lbf;
pub const Pressure = pressure_mod.Pressure;

pub const JPerKg = energy_per_mass_mod.JPerKg;
pub const BtuPerLbsm = energy_per_mass_mod.BtuPerLbsm;
pub const EnergyPerMass = energy_per_mass_mod.EnergyPerMass;

pub const JPerKgK = energy_per_mass_temperature_mod.JPerKgK;
pub const BtuPerLbsmR = energy_per_mass_temperature_mod.BtuPerLbsmR;
pub const EnergyPerMassTemperature = energy_per_mass_temperature_mod.EnergyPerMassTemperature;

pub const MPerSec = velocity_mod.MPerSec;
pub const FtPerSec = velocity_mod.FtPerSec;
pub const Velocity = velocity_mod.Velocity;

pub const M3PerKg = specific_volume_mod.M3PerKg;
pub const Ft3PerLbsm = specific_volume_mod.Ft3PerLbsm;
pub const SpecificVolume = specific_volume_mod.SpecificVolume;

pub const KgPerM3 = density_mod.KgPerM3;
pub const LbsmPerFt3 = density_mod.LbsmPerFt3;
pub const Density = density_mod.Density;

test "length conversion" {
    const one_meter = (Length{ .m = M.init(1.0) }).convertToSiUnit();
    const one_foot = (Length{ .ft = Ft.init(1.0) }).convertToSiUnit();
    const one_inch = (Length{ .inches = Inches.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, one_meter.value, 1e-9);
    try std.testing.expectApproxEqAbs(1.0 / feet_per_meter, one_foot.value, 1e-6);
    try std.testing.expectApproxEqAbs(1.0 / (feet_per_meter * inches_per_foot), one_inch.value, 1e-6);
}

test "area conversion" {
    const one_m2 = (Area{ .m2 = M2.init(1.0) }).convertToSiUnit();
    const one_ft2 = (Area{ .ft2 = Ft2.init(1.0) }).convertToSiUnit();
    const one_in2 = (Area{ .inches2 = Inches2.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, one_m2.value, 1e-9);
    try std.testing.expectApproxEqAbs(1.0 / feet_per_meter_sq, one_ft2.value, 1e-6);
    try std.testing.expectApproxEqAbs(1.0 / inches_per_meter_sq, one_in2.value, 1e-6);
}

test "volume conversion" {
    const one_m3 = (Volume{ .m3 = M3.init(1.0) }).convertToSiUnit();
    const one_ft3 = (Volume{ .ft3 = Ft3.init(1.0) }).convertToSiUnit();
    const one_in3 = (Volume{ .inches3 = Inches3.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, one_m3.value, 1e-9);
    try std.testing.expectApproxEqAbs(1.0 / feet_per_meter_cu, one_ft3.value, 1e-6);
    try std.testing.expectApproxEqAbs(1.0 / inches_per_meter_cu, one_in3.value, 1e-6);
}

test "volumetric flow rate conversion" {
    const per_sec = (VolumetricFlowRate{ .m3_per_sec = M3PerSec.init(1.0) }).convertToSiUnit();
    const per_min = (VolumetricFlowRate{ .m3_per_min = M3PerMin.init(1.0) }).convertToSiUnit();
    const ft_per_sec = (VolumetricFlowRate{ .ft3_per_sec = Ft3PerSec.init(1.0) }).convertToSiUnit();
    const ft_per_min = (VolumetricFlowRate{ .ft3_per_min = Ft3PerMin.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, per_sec.value, 1e-9);
    try std.testing.expectApproxEqAbs(60.0, per_min.value, 1e-6);
    try std.testing.expectApproxEqAbs(1.0 / feet_per_meter_cu, ft_per_sec.value, 1e-6);
    try std.testing.expectApproxEqAbs(60.0 / feet_per_meter_cu, ft_per_min.value, 1e-6);
}

test "mass conversion" {
    const kg_si = (Mass{ .kg = Kg.init(1.0) }).convertToSiUnit();
    const g_si = (Mass{ .g = G.init(1.0) }).convertToSiUnit();
    const lb_si = (Mass{ .lbsm = Lbsm.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, kg_si.value, 1e-9);
    try std.testing.expectApproxEqAbs(grams_per_kilogram, g_si.value, 1e-6);
    try std.testing.expectApproxEqAbs(pounds_mass_per_kilogram, lb_si.value, 1e-6);
}

test "temperature conversion" {
    const k_si = (Temperature{ .k = K.init(1.0) }).convertToSiUnit();
    const c_si = (Temperature{ .c = C.init(1.0) }).convertToSiUnit();
    const f_si = (Temperature{ .f = F.init(1.0) }).convertToSiUnit();
    const r_si = (Temperature{ .r = R.init(200.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, k_si.value, 1e-9);
    try std.testing.expectApproxEqAbs(274.15, c_si.value, 1e-6);
    try std.testing.expectApproxEqAbs(255.92777777777775, f_si.value, 1e-6);
    try std.testing.expectApproxEqAbs(111.11111111111111, r_si.value, 1e-6);
}

test "pressure conversion" {
    const pa_si = (Pressure{ .pa = Pa.init(1.0) }).convertToSiUnit();
    const kpa_si = (Pressure{ .kpa = KPa.init(1.0) }).convertToSiUnit();
    const lbf_si = (Pressure{ .lbf = Lbf.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, pa_si.value, 1e-9);
    try std.testing.expectApproxEqAbs(pascals_per_kpa, kpa_si.value, 1e-6);
    try std.testing.expectApproxEqAbs(pascals_per_psi, lbf_si.value, 1e-2);
}

test "energy per mass conversion" {
    const j_si = (EnergyPerMass{ .j_per_kg = JPerKg.init(1.0) }).convertToSiUnit();
    const btu_si = (EnergyPerMass{ .btu_per_lbsm = BtuPerLbsm.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, j_si.value, 1e-9);
    try std.testing.expectApproxEqAbs(energy_per_mass_conversion, btu_si.value, 1e-6);
}

test "energy per mass temperature conversion" {
    const j_per_si = (EnergyPerMassTemperature{ .j_per_kg_k = JPerKgK.init(1.0) }).convertToSiUnit();
    const btu_per_si = (EnergyPerMassTemperature{ .btu_per_lbsm_r = BtuPerLbsmR.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, j_per_si.value, 1e-9);
    try std.testing.expectApproxEqAbs(energy_per_mass_temp_conversion, btu_per_si.value, 1e-6);
}

test "velocity conversion" {
    const m_si = (Velocity{ .m_per_sec = MPerSec.init(1.0) }).convertToSiUnit();
    const ft_si = (Velocity{ .ft_per_sec = FtPerSec.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, m_si.value, 1e-9);
    try std.testing.expectApproxEqAbs(1.0 / feet_per_meter, ft_si.value, 1e-6);
}

test "specific volume conversion" {
    const si = (SpecificVolume{ .m3_per_kg = M3PerKg.init(1.0) }).convertToSiUnit();
    const imperial = (SpecificVolume{ .ft3_per_lbsm = Ft3PerLbsm.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, si.value, 1e-9);
    try std.testing.expectApproxEqAbs(specific_volume_conversion, imperial.value, 1e-6);
}

test "density conversion" {
    const si = (Density{ .kg_per_m3 = KgPerM3.init(1.0) }).convertToSiUnit();
    const imperial = (Density{ .lbsm_per_ft3 = LbsmPerFt3.init(1.0) }).convertToSiUnit();

    try std.testing.expectApproxEqAbs(1.0, si.value, 1e-9);
    try std.testing.expectApproxEqAbs(density_conversion, imperial.value, 1e-6);
}
