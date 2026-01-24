const common = @import("common.zig");

const UnitValueType = common.UnitValueType;
const UnitLabel = common.UnitLabel;
const RawUnit = common.RawUnit;
const ParseUnitError = common.ParseUnitError;
const isMatch = common.isMatch;
const watts_per_kw = common.watts_per_kw;
const watts_per_hp = common.watts_per_hp;
const watts_per_btu_hr = common.watts_per_btu_hr;

pub const W = UnitValueType();
pub const kW = UnitValueType();
pub const Hp = UnitValueType();
pub const BtuPerHr = UnitValueType();

const power_labels = [_]UnitLabel{
    .{ .abbreviation = "W", .plural = "watts" },
    .{ .abbreviation = "kW", .plural = "kilowatts" },
    .{ .abbreviation = "hp", .plural = "horsepower" },
    .{ .abbreviation = "BTU/hr", .plural = "british thermal units per hour" },
};

pub const Power = union(enum) {
    w: W,
    kw: kW,
    hp: Hp,
    btu_per_hr: BtuPerHr,

    pub fn convertToSiUnit(self: Power) W {
        return switch (self) {
            .w => |val| val,
            .kw => |val| W.init(val.value * watts_per_kw),
            .hp => |val| W.init(val.value * watts_per_hp),
            .btu_per_hr => |val| W.init(val.value * watts_per_btu_hr),
        };
    }

    pub fn listUnitLabels() []const UnitLabel {
        return &power_labels;
    }

    pub fn getSiUnitLabel() UnitLabel {
        return power_labels[0];
    }

    pub fn getValue(self: Power) f64 {
        return switch (self) {
            .w => |val| val.value,
            .kw => |val| val.value,
            .hp => |val| val.value,
            .btu_per_hr => |val| val.value,
        };
    }

    pub fn tryConvert(self: Power, unit_display: []const u8) ParseUnitError!Power {
        const value_si = self.convertToSiUnit().value;

        if (isMatch(unit_display, "W")) return .{ .w = W.init(value_si) };
        if (isMatch(unit_display, "kW")) return .{ .kw = kW.init(value_si / watts_per_kw) };
        if (isMatch(unit_display, "hp")) return .{ .hp = Hp.init(value_si / watts_per_hp) };
        if (isMatch(unit_display, "BTU/hr")) return .{ .btu_per_hr = BtuPerHr.init(value_si / watts_per_btu_hr) };

        return ParseUnitError.UnknownUnit;
    }

    pub fn fromRawUnit(raw: RawUnit) ParseUnitError!Power {
        if (isMatch(raw.unit_display, "W")) return .{ .w = W.init(raw.value) };
        if (isMatch(raw.unit_display, "kW")) return .{ .kw = kW.init(raw.value) };
        if (isMatch(raw.unit_display, "hp")) return .{ .hp = Hp.init(raw.value) };
        if (isMatch(raw.unit_display, "BTU/hr")) return .{ .btu_per_hr = BtuPerHr.init(raw.value) };
        return ParseUnitError.UnknownUnit;
    }

    pub fn toRawUnit(self: Power) RawUnit {
        return switch (self) {
            .w => |val| .{ .value = val.value, .unit_display = "W" },
            .kw => |val| .{ .value = val.value, .unit_display = "kW" },
            .hp => |val| .{ .value = val.value, .unit_display = "hp" },
            .btu_per_hr => |val| .{ .value = val.value, .unit_display = "BTU/hr" },
        };
    }
};
