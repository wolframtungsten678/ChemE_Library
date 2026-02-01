const std = @import("std");
const iapws97 = @import("iapws97.zig");
const units = @import("../../units.zig");

pub const RankineCycleArgs = struct {
    boilerTemperature: units.Temperature,
    boilerPressure: units.Pressure,
    condenserPressure: units.Pressure,
    pumpEfficiency: f64,
    turbineEfficiency: f64,
    powerRequirement: units.Power,
};

pub const RankineCycleResult = struct {
    condenserSteamQuality: f64,
    boilerWork: units.EnergyPerMass,
    turbineWork: units.EnergyPerMass,
    thermalEfficiency: f64,
    pumpWork: units.EnergyPerMass,
    condenserWork: units.EnergyPerMass,
    netWork: units.EnergyPerMass,
    steamRate: units.MassFlowRate,
    boilerHeatTransferRate: units.MassFlowRate,
    condenserHeatTransferRate: units.MassFlowRate,
};

// force error, but bypass compiletime check
fn kill() !void {
    var x: u64 = 0;
    if (x > 1) {
        x += 1;
    } else {
        return error.Kill;
    }
}

pub fn rankineCycle(args: RankineCycleArgs) !RankineCycleResult {
    const boiler_conditions = try iapws97.getSteamEntryByPressureAndTemperature(args.boilerPressure, args.boilerTemperature);
    std.debug.print("{any}\n", .{boiler_conditions.temperature});
    std.debug.print("{any}\n", .{boiler_conditions.pressure});
    std.debug.print("{any}\n", .{boiler_conditions.enthalpy});
    std.debug.print("{any}\n", .{boiler_conditions.entropy});
    const condenser_liquid_conditions = try iapws97.getSteamEntryBySatPressure(iapws97.SteamNonCriticalPhaseRegion.Liquid, args.condenserPressure);
    std.debug.print("{any}\n", .{condenser_liquid_conditions.enthalpy});

    const ideal_condenser_conditions = try iapws97.getSteamEntryByPressureAndEntropy(condenser_liquid_conditions.pressure, boiler_conditions.entropy);
    std.debug.print("{any}\n", .{ideal_condenser_conditions.enthalpy});

    const boiler_enthalpy_raw = boiler_conditions.enthalpy.convertToSiUnit().value;
    const ideal_condenser_enthalpy_raw = ideal_condenser_conditions.enthalpy.convertToSiUnit().value;
    const ideal_energy_released = boiler_enthalpy_raw - ideal_condenser_enthalpy_raw;
    const actual_energy_released = args.turbineEfficiency * ideal_energy_released;
    std.debug.print("{any}\n", .{actual_energy_released});

    const condenser_pressure_raw = condenser_liquid_conditions.pressure.convertToSiUnit().value;
    const actual_condenser_enthalpy_raw = boiler_enthalpy_raw - actual_energy_released;
    const actual_condenser_enthalpy = units.EnergyPerMass{ .j_per_kg = units.JPerKg.init(actual_condenser_enthalpy_raw) };
    const actual_condenser_conditions = try iapws97.getSteamEntryByPressureAndEnthalpy(condenser_liquid_conditions.pressure, actual_condenser_enthalpy);

    const turbine_work_raw = boiler_enthalpy_raw - actual_condenser_enthalpy_raw;
    const turbine_work = units.EnergyPerMass{ .j_per_kg = units.JPerKg.init(turbine_work_raw) };
    //const net_work_raw = turbine_work_raw - pump_work_raw;
    const net_work_raw = turbine_work_raw;
    const net_work = units.EnergyPerMass{ .j_per_kg = units.JPerKg.init(net_work_raw) };
    //const thermal_efficiency = @abs(turbine_work_raw) / boiler_work_raw;
    const thermal_efficiency: i64 = 0;
    const power_requirement_raw = args.powerRequirement.convertToSiUnit().value;

    const steam_rate_raw = power_requirement_raw / net_work_raw;
    const steam_rate = units.MassFlowRate{ .kg_per_sec = units.KgPerSec.init(steam_rate_raw) };
    std.debug.print("{any}\n", .{actual_condenser_conditions.phase_region});
    std.debug.print("{any}\n", .{actual_condenser_conditions.entropy});
    std.debug.print("{any}\n", .{args.powerRequirement});
    std.debug.print("{any}\n", .{steam_rate});

    try kill();

    const ideal_condenser_steam_quality = switch (ideal_condenser_conditions.phase_region) {
        .Composite => |composite| composite.liquid_vapor.vapor_frac,
        // TODO: make this a more appropriate error
        else => return iapws97.SteamError.InvalidPhaseFractions,
    };
    std.debug.print("{any}\n", .{ideal_condenser_steam_quality});

    const condenser_sv = condenser_liquid_conditions.specific_volume.convertToSiUnit().value;
    const boiler_pressure = boiler_conditions.pressure.convertToSiUnit().value;
    const pump_efficiency = args.pumpEfficiency;
    if (pump_efficiency < 0.0 or pump_efficiency > 1.0) return error.InvalidPumpEfficiency;

    const pump_work_raw = ((condenser_sv * (boiler_pressure - condenser_pressure_raw))) / pump_efficiency;

    const pump_work = units.EnergyPerMass{ .j_per_kg = units.JPerKg.init(pump_work_raw) };
    const pump_enthalpy_raw = ideal_condenser_enthalpy_raw + pump_work_raw;
    std.debug.print("pump enthalpy:\n{}\n", .{pump_enthalpy_raw});

    const boiler_work_raw = boiler_enthalpy_raw - pump_enthalpy_raw;
    std.debug.print("Boiler Work:\n{} = {} - {}\n", .{ boiler_work_raw, boiler_enthalpy_raw, pump_enthalpy_raw });

    const boiler_work = units.EnergyPerMass{ .j_per_kg = units.JPerKg.init(boiler_work_raw) };

    const turbine_efficiency = args.turbineEfficiency;
    if (turbine_efficiency < 0.0 or turbine_efficiency > 1.0) return error.InvalidTurbineEfficiency;

    //const condenser_work_raw = -(condenser_enthalpy_raw - (boiler_enthalpy_raw - turbine_work_raw));
    const condenser_work_raw = -(ideal_condenser_enthalpy_raw - (boiler_enthalpy_raw - turbine_work_raw));
    std.debug.print("Steam Quality:\n{}\n{any}\n", .{ ideal_condenser_steam_quality, ideal_condenser_conditions.enthalpy });
    std.debug.print("Condenser Work:\n{} = -({} - ({} - {}))\n", .{ condenser_work_raw, ideal_condenser_enthalpy_raw, boiler_enthalpy_raw, turbine_work_raw });
    std.debug.print("Condenser Work:\n{} = -({} - {})\n", .{ condenser_work_raw, ideal_condenser_enthalpy_raw, boiler_enthalpy_raw - turbine_work_raw });
    const condenser_work = units.EnergyPerMass{ .j_per_kg = units.JPerKg.init(condenser_work_raw) };

    const boiler_heat_transfer_rate_raw = steam_rate_raw * boiler_work_raw;
    const boiler_heat_transfer_rate = units.MassFlowRate{ .kg_per_sec = units.KgPerSec.init(boiler_heat_transfer_rate_raw) };

    const condenser_heat_transfer_rate_raw = steam_rate_raw * condenser_work_raw;
    const condenser_heat_transfer_rate = units.MassFlowRate{ .kg_per_sec = units.KgPerSec.init(condenser_heat_transfer_rate_raw) };

    return RankineCycleResult{
        .condenserSteamQuality = ideal_condenser_steam_quality,
        .boilerWork = boiler_work,
        .turbineWork = turbine_work,
        .thermalEfficiency = thermal_efficiency,
        .netWork = net_work,
        .condenserWork = condenser_work,
        .pumpWork = pump_work,
        .steamRate = steam_rate,
        .boilerHeatTransferRate = boiler_heat_transfer_rate,
        .condenserHeatTransferRate = condenser_heat_transfer_rate,
    };
}

// TODO: Add a test case for this:
// BoilerTemperature.Value = 500;
// BoilerPressure.Value = 8600e3;
// CondenserPressure.Value = 10e3;
// PumpEfficiency.Value = 0.75;
// TurbineEfficiency.Value = 0.75;
// PowerRequirement.Value = 80e3;
test "rankine cycle example inputs" {
    const args = RankineCycleArgs{
        // TODO: Something like this should not compile
        // .c and K should not work together
        // .boilerTemperature = units.Temperature{ .c = units.K.init(500.0) },
        .boilerTemperature = units.Temperature{ .c = units.C.init(500.0) },
        .boilerPressure = units.Pressure{ .pa = units.Pa.init(8600e3) },
        .condenserPressure = units.Pressure{ .pa = units.Pa.init(10e3) },
        .pumpEfficiency = 0.75,
        .turbineEfficiency = 0.75,
        .powerRequirement = units.Power{ .kw = units.kW.init(80e3) },
    };

    const result = try rankineCycle(args);

    const condenser_steam_quality = result.condenserSteamQuality;
    const pump_work = result.pumpWork.convertToSiUnit().value;
    const boiler_work = result.boilerWork.convertToSiUnit().value;
    const condenser_work = result.condenserWork.convertToSiUnit().value;
    const turbine_work = result.turbineWork.convertToSiUnit().value;
    const thermal_efficiency = result.thermalEfficiency;
    const net_work = result.netWork.convertToSiUnit().value;
    const steam_rate = result.steamRate.convertToSiUnit().value;
    const boiler_heat_rate = result.boilerHeatTransferRate.convertToSiUnit().value;
    const condenser_heat_rate = result.condenserHeatTransferRate.convertToSiUnit().value;

    try std.testing.expectApproxEqAbs(0.8051, condenser_steam_quality, 1e-3);
    try std.testing.expectApproxEqAbs(11570, pump_work, 1);
    try std.testing.expectApproxEqAbs(3189e3, boiler_work, 1e3);
    try std.testing.expectApproxEqAbs(0, turbine_work, 1e-3);
    // check this again
    try std.testing.expectApproxEqAbs(0, condenser_work, 1);
    try std.testing.expectApproxEqAbs(0, turbine_work, 1e-3);
    try std.testing.expectApproxEqAbs(0, thermal_efficiency, 1e-3);
    try std.testing.expectApproxEqAbs(0, net_work, 1e-3);
    try std.testing.expectApproxEqAbs(0, steam_rate, 1e-3);
    try std.testing.expectApproxEqAbs(0, boiler_heat_rate, 1e-3);
    try std.testing.expectApproxEqAbs(0, condenser_heat_rate, 1e-3);
}
