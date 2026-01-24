const iapws97 = @import("iapws97.zig");
const units = @import("../../units.zig");

pub const RankineCycleArgs = struct {
    boilerTemperature: units.Temperature,
    boilerPressure: units.Pressure,
    condenserPressure: units.Pressure,
    pump_efficiency: f64,
};

pub const RankineCycleResult = struct {
    condenserSteamQuality: f64,
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
    const pump_efficiency = args.pump_efficiency;
    if (pump_efficiency < 0.0 or pump_efficiency > 1.0) return error.InvalidPumpEfficiency;

    const pump_work_raw = ((condenser_sv * (boiler_pressure - condenser_pressure)) * 1e-3) / pump_efficiency;

    const pump_work = units.EnergyPerMass(units.JPerKg.init(pump_work_raw));
    const condenser_enthalpy_raw = condenser_liquid_conditions.enthalpy.convertToSiUnit().value;
    const boiler_enthalpy_raw = condenser_enthalpy_raw + pump_work_raw;

    const boiler_work = boiler_enthalpy_raw - boiler_enthalpy_raw;

    //// in kj / kg
    //double boilerEnthalpy = condenserLiquidConditions.H + PumpWork.Value;

    //BoilerWork.Value = boilerConditions.H - boilerEnthalpy;

    //TurbineWork.Value = -(condenserConditions.H - boilerConditions.H) * TurbineEfficiency.Value;

    //CondenserWork.Value = -(condenserLiquidConditions.H - (boilerConditions.H - TurbineWork.Value));

    //NetWork.Value =  TurbineWork.Value - PumpWork.Value;

    //ThermalEfficiency.Value = Math.Abs(TurbineWork.Value) / BoilerWork.Value;

    //SteamRate.Value = PowerRequirement.Value / NetWork.Value;

    //BoilerHeatTransRate.Value = SteamRate.Value * BoilerWork.Value;

    //CondenserHeatTransRate.Value = SteamRate.Value * CondenserWork.Value;

    return .{
        .condenserSteamQuality = condenserSteamQuality,
    };
}
