#include "Air_control.h"
#include <cmath>
#include <iostream>

namespace AirControl {
    // Константы из ГОСТ 4401-81
    constexpr double R_EARTH = 6356767.0;        // Условный радиус Земли, м
    constexpr double G_SEA_LEVEL = 9.80665;      // Ускорение свободного падения на уровне моря, м/с²
    constexpr double P_SEA_LEVEL = 101325.0;     // Давление на уровне моря, Па
    constexpr double T_SEA_LEVEL = 288.15;       // Температура на уровне моря, K
    constexpr double R_AIR = 287.05287;          // Газовая постоянная для воздуха, Дж/(кг·K)
    constexpr double GAMMA = 1.4;                // Показатель адиабаты для воздуха

    double AtmosphereCalculator::GeometricToGeopotential(double h) {
        // Формула пересчета геометрической высоты в геопотенциальную
        // H = (R * h) / (R + h)
        return (R_EARTH * h) / (R_EARTH + h);
    }

    double AtmosphereCalculator::CalculateGravity(double geometric_altitude) {
        // Ускорение свободного падения зависит от ГЕОМЕТРИЧЕСКОЙ высоты
        double ratio = R_EARTH / (R_EARTH + geometric_altitude);
        return G_SEA_LEVEL * ratio * ratio;
    }

    double AtmosphereCalculator::CalculateSpeedOfSound(double T) {
        return std::sqrt(GAMMA * R_AIR * T);
    }

    // ВСЕ расчеты параметров атмосферы должны использовать ГЕОПОТЕНЦИАЛЬНУЮ высоту
    double AtmosphereCalculator::CalculateTemperature(double H_geopotential) {
        // H_geopotential - геопотенциальная высота в метрах

        if (H_geopotential < -2000) {
            return T_SEA_LEVEL - 0.0065 * H_geopotential;
        }
        else if (H_geopotential <= 0) {
            return T_SEA_LEVEL - 0.0065 * H_geopotential;
        }
        else if (H_geopotential <= 11000) {
            return T_SEA_LEVEL - 0.0065 * H_geopotential;
        }
        else if (H_geopotential <= 20000) {
            return 216.65;
        }
        else if (H_geopotential <= 32000) {
            return 216.65 + 0.0010 * (H_geopotential - 20000);
        }
        else if (H_geopotential <= 47000) {
            return 228.65 + 0.0028 * (H_geopotential - 32000);
        }
        else if (H_geopotential <= 51000) {
            return 270.65;
        }
        else if (H_geopotential <= 71000) {
            return 270.65 - 0.0028 * (H_geopotential - 51000);
        }
        else if (H_geopotential <= 85000) {
            return 214.65 - 0.0020 * (H_geopotential - 71000);
        }
        else {
            return 186.65;
        }
    }

    double AtmosphereCalculator::CalculatePressure(double H_geopotential, double T) {
        // Расчет давления использует ГЕОПОТЕНЦИАЛЬНУЮ высоту

        if (H_geopotential <= 0) {
            double beta = -0.0065;
            return P_SEA_LEVEL * std::pow(T_SEA_LEVEL / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H_geopotential <= 11000) {
            double beta = -0.0065;
            return P_SEA_LEVEL * std::pow(T_SEA_LEVEL / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H_geopotential <= 20000) {
            double P_11000 = CalculatePressure(11000, 216.65);
            return P_11000 * std::exp(-G_SEA_LEVEL * (H_geopotential - 11000) / (R_AIR * 216.65));
        }
        else if (H_geopotential <= 32000) {
            double P_20000 = CalculatePressure(20000, 216.65);
            double beta = 0.0010;
            double T_20000 = 216.65;
            return P_20000 * std::pow(T_20000 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H_geopotential <= 47000) {
            double P_32000 = CalculatePressure(32000, 228.65);
            double beta = 0.0028;
            double T_32000 = 228.65;
            return P_32000 * std::pow(T_32000 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H_geopotential <= 51000) {
            double P_47000 = CalculatePressure(47000, 270.65);
            return P_47000 * std::exp(-G_SEA_LEVEL * (H_geopotential - 47000) / (R_AIR * 270.65));
        }
        else if (H_geopotential <= 71000) {
            double P_51000 = CalculatePressure(51000, 270.65);
            double beta = -0.0028;
            double T_51000 = 270.65;
            return P_51000 * std::pow(T_51000 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H_geopotential <= 85000) {
            double P_71000 = CalculatePressure(71000, 214.65);
            double beta = -0.0020;
            double T_71000 = 214.65;
            return P_71000 * std::pow(T_71000 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else {
            double P_85000 = CalculatePressure(85000, 186.65);
            return P_85000 * std::exp(-G_SEA_LEVEL * (H_geopotential - 85000) / (R_AIR * 186.65));
        }
    }

    double AtmosphereCalculator::CalculateDensity(double P, double T) {
        return P / (R_AIR * T);
    }

    AtmosphereParams AtmosphereCalculator::Calculate(double geometric_altitude) {
        AtmosphereParams params;

        // Сохраняем исходную геометрическую высоту
        params.geometric_altitude = geometric_altitude;

        // Преобразуем в геопотенциальную высоту
        params.geopotential_altitude = GeometricToGeopotential(geometric_altitude);

        // ВСЕ расчеты параметров атмосферы используют геопотенциальную высоту
        params.temperature = CalculateTemperature(params.geopotential_altitude);
        params.pressure = CalculatePressure(params.geopotential_altitude, params.temperature);
        params.density = CalculateDensity(params.pressure, params.temperature);
        params.speed_of_sound = CalculateSpeedOfSound(params.temperature);

        // Только ускорение свободного падения использует геометрическую высоту
        params.gravity = CalculateGravity(geometric_altitude);

        return params;
    }

    std::vector<AtmosphereParams> AtmosphereCalculator::CalculateRange(double min_altitude, double max_altitude, double step) {
        std::vector<AtmosphereParams> results;

        for (double h = min_altitude; h <= max_altitude; h += step) {
            results.push_back(Calculate(h));
        }

        return results;
    }
}