#pragma once
#include <vector>
#include <string>

namespace AirControl {
    struct AtmosphereParams {
        double geometric_altitude;    // Геометрическая высота, м
        double geopotential_altitude; // Геопотенциальная высота, м
        double temperature;           // Температура, K
        double pressure;              // Давление, Па
        double density;               // Плотность, кг/м³
        double speed_of_sound;        // Скорость звука, м/с
        double gravity;               // Ускорение свободного падения, м/с²
    };

    class AtmosphereCalculator {
    public:
        // Расчет параметров для диапазона высот
        static std::vector<AtmosphereParams> CalculateRange(double min_altitude, double max_altitude, double step = 100.0);

        // Расчет параметров для одной высоты
        static AtmosphereParams Calculate(double geometric_altitude);

    private:
        // Вспомогательные функции
        static double GeometricToGeopotential(double h);
        static double CalculateGravity(double geometric_altitude);
        static double CalculateTemperature(double H);
        static double CalculatePressure(double H, double T);
        static double CalculateDensity(double P, double T);
        static double CalculateSpeedOfSound(double T);  // Новая функция для расчета скорости звука
    };
}