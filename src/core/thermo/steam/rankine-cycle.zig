const iapws97 = @import("iapws97.zig");
const units = @import("../../units.zig");

pub const RankineCycleArgs = struct {
    // BoilerTemperature.Value = 500;
    // BoilerPressure.Value = 8600e3;
    // CondenserPressure.Value = 10e3;
    // PumpEfficiency.Value = 0.75;
    // TurbineEfficiency.Value = 0.75;
    // PowerRequirement.Value = 80e3;
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
};

pub fn rankineCycle(args: RankineCycleArgs) !RankineCycleResult {
    const boiler_conditions = try iapws97.getSteamEntryByPressureAndTemperature(args.boilerPressure, args.boilerTemperature);
    const condenser_liquid_conditions = try iapws97.getSteamEntryBySatPressure(args.condenserPressure, iapws97.NonCriticalPhaseRegion.Liquid);

    const condenser_conditions = try iapws97.getSteamEntryByPressureAndEntropy(boiler_conditions.pressure, boiler_conditions.entropy);

    const condenser_steam_quality = switch (condenser_conditions.phase_region) {
        .Composite => |composite| composite.liquid_vapor.vapor_frac,
        else => return iapws97.SteamError.InvalidPhaseFractions,
    };

    const condenser_sv = condenser_liquid_conditions.specific_volume.convertToSiUnit().value;
    const boiler_pressure = boiler_conditions.pressure.convertToSiUnit().value;
    const condenser_pressure = condenser_liquid_conditions.pressure.convertToSiUnit().value;
    const pump_efficiency = args.pumpEfficiency;
    if (pump_efficiency < 0.0 or pump_efficiency > 1.0) return error.InvalidPumpEfficiency;

    const pump_work_raw = ((condenser_sv * (boiler_pressure - condenser_pressure)) * 1e-3) / pump_efficiency;

    const pump_work = units.EnergyPerMass(units.JPerKg.init(pump_work_raw));
    const condenser_enthalpy_raw = condenser_liquid_conditions.enthalpy.convertToSiUnit().value;
    // TODO: make this a better name
    const boiler_enthalpy_raw_better_name = condenser_enthalpy_raw + pump_work_raw;

    const boiler_enthalpy_raw = boiler_conditions.enthalpy.convertToSiUnit().value;
    const boiler_work_raw = boiler_enthalpy_raw - boiler_enthalpy_raw_better_name;

    const boiler_work = units.EnergyPerMass(units.JPerKg.init(boiler_work_raw));

    const turbine_efficiency = args.turbineEfficiency;
    if (turbine_efficiency < 0.0 or turbine_efficiency > 1.0) return error.InvalidTurbineEfficiency;

    const turbine_work_raw = -(condenser_enthalpy_raw - boiler_enthalpy_raw) * turbine_efficiency;
    const turbine_work = units.EnergyPerMass(units.JPerKg.init(turbine_work_raw));

    const condenser_work_raw = -(condenser_enthalpy_raw - (boiler_enthalpy_raw - turbine_work_raw));
    const condenser_work = units.EnergyPerMass(units.JPerKg.init(condenser_work_raw));

    const net_work_raw = turbine_work_raw - pump_work_raw;
    const net_work = units.EnergyPerMass(units.JPerKg.init(net_work_raw));
    const thermal_efficiency = @abs(turbine_work_raw) / boiler_work_raw;
    const power_requirement_raw = args.powerRequirement.convertToSiUnit().value;
    // TODO: make a mass flow rate unit
    const steam_rate_raw = power_requirement_raw / net_work_raw;

    const boiler_heat_transfer_rate_raw = steam_rate_raw * boiler_work_raw;
    _ = boiler_heat_transfer_rate_raw;

    const condenser_heat_transfer_rate_raw = steam_rate_raw * condenser_work_raw;
    _ = condenser_heat_transfer_rate_raw;

    return RankineCycleResult{
        .condenserSteamQuality = condenser_steam_quality,
        .boilerWork = boiler_work,
        .turbineWork = turbine_work,
        .thermalEfficiency = thermal_efficiency,
        .netWork = net_work,
        .condenserWork = condenser_work,
        .pumpWork = pump_work,
    };
}
