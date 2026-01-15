#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <clocale> // Для setlocale
#include "../Air_control/Air_control.h"

using namespace std;

// Определяем константу PI, если она не определена
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Функция для настройки русской кодировки консоли
void SetConsoleToRussian() {
    setlocale(LC_ALL, "Russian");
    // Для консоли Windows устанавливаем кодовую страницу 866 (DOS) или 1251 (Windows)
    system("chcp 1251 > nul"); // Устанавливаем кодовую страницу Windows-1251
}

// Константы и параметры ЛА
const double m0 = 1255.0;        // начальная масса, кг
const double Jz = 127.0;         // момент инерции, кг*м^2
const double l_d = 0.215;        // запас устойчивости, м
const double S_a = 0.14;         // площадь выходного сечения сопла, м^2
const double S_m = 0.231;        // характерная площадь, м^2
const double g0 = 9.80665;       // ускорение свободного падения на уровне моря, м/с^2
const double p0N = 101325.0;     // нормальное давление у поверхности Земли, Па

// Исходные данные для варианта Трояновский (15)
const double V0 = 65.0;          // начальная скорость, м/с
const double theta_c0_deg = 48.0; // начальный угол наклона траектории, градусы
const double theta_c0 = theta_c0_deg * M_PI / 180.0; // в радианах
const double m_dot = 85.0;       // массовый расход, кг/с
const double W = 2180.0;         // скорость истечения газов, м/с
const double y0_val = 3869.0;    // начальная высота, м (переименовано, чтобы избежать конфликта)
const double omega_z0 = 0.03;    // начальная угловая скорость, с^-1
const double theta0 = theta_c0;  // начальный угол тангажа, рад
const double t0 = 0.0;           // начальное время, с
const double t_k = 3.72;         // время активного участка, с
const double dt = 0.1;           // шаг интегрирования, с

// Табличные данные для Cxa(M) и Cya(M)
vector<double> M_values = { 0.01, 0.55, 0.8, 0.9, 1.0, 1.06, 1.1, 1.2, 1.3, 1.4, 2.0, 2.6, 3.4, 6.0, 10.2 };
vector<double> Cxa_values = { 0.30, 0.30, 0.55, 0.70, 0.84, 0.86, 0.87, 0.83, 0.80, 0.79, 0.65, 0.55, 0.50, 0.45, 0.41 };
vector<double> Cya_values = { 0.25, 0.25, 0.25, 0.20, 0.30, 0.31, 0.25, 0.25, 0.25, 0.25, 0.20, 0.25, 0.25, 0.25, 0.15 };

// Функция линейной интерполяции
double interpolate(const vector<double>& x_vals, const vector<double>& y_vals, double x) {
    if (x <= x_vals[0]) return y_vals[0];
    if (x >= x_vals.back()) return y_vals.back();
    for (size_t i = 0; i < x_vals.size() - 1; ++i) {
        if (x >= x_vals[i] && x <= x_vals[i + 1]) {
            return y_vals[i] + (y_vals[i + 1] - y_vals[i]) * (x - x_vals[i]) / (x_vals[i + 1] - x_vals[i]);
        }
    }
    return y_vals.back();
}

// Функция для расчёта Cxa(M)
double get_Cxa(double M) {
    return interpolate(M_values, Cxa_values, M);
}

// Функция для расчёта Cya(M)
double get_Cya(double M) {
    return interpolate(M_values, Cya_values, M);
}

// Структура для хранения состояния системы
struct State {
    double V;      // скорость, м/с
    double theta_c; // угол наклона траектории, рад
    double x_c;    // координата x, м
    double y_c;    // координата y, м
    double omega_z; // угловая скорость, рад/с
    double theta;  // угол тангажа, рад
    double t;      // время, с
    double m;      // масса, кг
};

// Функция для вычисления производных (правые части уравнений)
void derivatives(const State& s, State& dsdt, const AirControl::AtmosphereParams& atm) {
    double alpha = s.theta - s.theta_c; // угол атаки, рад
    double M = s.V / atm.speed_of_sound; // число Маха

    // Аэродинамические коэффициенты
    double Cxa = get_Cxa(M);
    double Cya_alpha = get_Cya(M);
    double Cya = Cya_alpha * alpha;

    // Динамическое давление
    double q = 0.5 * atm.density * s.V * s.V;

    // Тяга
    double P0 = m_dot * W;
    double P = P0 + S_a * (p0N - atm.pressure);

    // Аэродинамические силы
    double X_a = Cxa * S_m * q;
    double Y_a = Cya * S_m * q;

    // Проекции веса
    double Q_xk = s.m * atm.gravity * sin(s.theta_c);
    double Q_yk = s.m * atm.gravity * cos(s.theta_c);

    // Суммарные силы
    double F_xk = P * cos(alpha) - X_a - Q_xk;
    double F_yk = P * sin(alpha) + Y_a - Q_yk;

    // Момент
    double M_z = -(Y_a * cos(alpha) + X_a * sin(alpha)) * l_d;

    // Производные
    dsdt.V = F_xk / s.m;
    if (abs(s.V) > 1e-6) {
        dsdt.theta_c = F_yk / (s.m * s.V);
    }
    else {
        dsdt.theta_c = 0.0;
    }
    dsdt.x_c = s.V * cos(s.theta_c);
    dsdt.y_c = s.V * sin(s.theta_c);
    dsdt.omega_z = M_z / Jz;
    dsdt.theta = s.omega_z;
    dsdt.t = 1.0;
    dsdt.m = -m_dot;
}

// Метод Рунге-Кутта 4-го порядка
void rungeKuttaStep(State& s, double dt, const AirControl::AtmosphereParams& atm) {
    State k1, k2, k3, k4, temp;

    // k1
    derivatives(s, k1, atm);

    // k2
    temp.V = s.V + 0.5 * dt * k1.V;
    temp.theta_c = s.theta_c + 0.5 * dt * k1.theta_c;
    temp.x_c = s.x_c + 0.5 * dt * k1.x_c;
    temp.y_c = s.y_c + 0.5 * dt * k1.y_c;
    temp.omega_z = s.omega_z + 0.5 * dt * k1.omega_z;
    temp.theta = s.theta + 0.5 * dt * k1.theta;
    temp.t = s.t + 0.5 * dt;
    temp.m = s.m + 0.5 * dt * k1.m;
    derivatives(temp, k2, atm);

    // k3
    temp.V = s.V + 0.5 * dt * k2.V;
    temp.theta_c = s.theta_c + 0.5 * dt * k2.theta_c;
    temp.x_c = s.x_c + 0.5 * dt * k2.x_c;
    temp.y_c = s.y_c + 0.5 * dt * k2.y_c;
    temp.omega_z = s.omega_z + 0.5 * dt * k2.omega_z;
    temp.theta = s.theta + 0.5 * dt * k2.theta;
    temp.t = s.t + 0.5 * dt;
    temp.m = s.m + 0.5 * dt * k2.m;
    derivatives(temp, k3, atm);

    // k4
    temp.V = s.V + dt * k3.V;
    temp.theta_c = s.theta_c + dt * k3.theta_c;
    temp.x_c = s.x_c + dt * k3.x_c;
    temp.y_c = s.y_c + dt * k3.y_c;
    temp.omega_z = s.omega_z + dt * k3.omega_z;
    temp.theta = s.theta + dt * k3.theta;
    temp.t = s.t + dt;
    temp.m = s.m + dt * k3.m;
    derivatives(temp, k4, atm);

    // Обновление состояния
    s.V += dt / 6.0 * (k1.V + 2 * k2.V + 2 * k3.V + k4.V);
    s.theta_c += dt / 6.0 * (k1.theta_c + 2 * k2.theta_c + 2 * k3.theta_c + k4.theta_c);
    s.x_c += dt / 6.0 * (k1.x_c + 2 * k2.x_c + 2 * k3.x_c + k4.x_c);
    s.y_c += dt / 6.0 * (k1.y_c + 2 * k2.y_c + 2 * k3.y_c + k4.y_c);
    s.omega_z += dt / 6.0 * (k1.omega_z + 2 * k2.omega_z + 2 * k3.omega_z + k4.omega_z);
    s.theta += dt / 6.0 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
    s.t += dt;
    s.m += dt / 6.0 * (k1.m + 2 * k2.m + 2 * k3.m + k4.m);

    // Защита от отрицательной массы
    if (s.m < 0) s.m = 0;
}

int main() {
    // Настройка русской кодировки консоли
    SetConsoleToRussian();

    cout << "==========================================" << endl;
    cout << "Расчет траектории ЛА (вариант Трояновский)" << endl;
    cout << "==========================================" << endl;
    cout << "Исходные данные:" << endl;
    cout << "V0 = " << V0 << " м/с" << endl;
    cout << "theta_c0 = " << theta_c0_deg << " град" << endl;
    cout << "m_dot = " << m_dot << " кг/с" << endl;
    cout << "W = " << W << " м/с" << endl;
    cout << "y0 = " << y0_val << " м" << endl;
    cout << "omega_z0 = " << omega_z0 << " 1/с" << endl;
    cout << "t_k = " << t_k << " с" << endl;
    cout << "dt = " << dt << " с" << endl;
    cout << "==========================================" << endl;

    // Начальное состояние
    State s;
    s.V = V0;
    s.theta_c = theta_c0;
    s.x_c = 0.0;
    s.y_c = y0_val;
    s.omega_z = omega_z0;
    s.theta = theta0;
    s.t = t0;
    s.m = m0;

    // Вектор для хранения результатов
    vector<State> results;
    results.push_back(s);

    // Интегрирование по времени
    int step_count = 0;
    cout << "\nНачало интегрирования..." << endl;

    for (double t = t0 + dt; t <= t_k + 0.5 * dt; t += dt) {
        step_count++;

        // Расчет параметров атмосферы на текущей высоте
        AirControl::AtmosphereParams atm = AirControl::AtmosphereCalculator::Calculate(s.y_c);

        // Шаг метода Рунге-Кутта
        rungeKuttaStep(s, dt, atm);

        // Сохранение результатов
        results.push_back(s);

        // Вывод прогресса
        if (step_count % 10 == 0) {
            cout << "Шаг " << step_count << ": t = " << s.t << " с, y = " << s.y_c << " м, V = " << s.V << " м/с" << endl;
        }
    }

    cout << "\nИнтегрирование завершено. Всего шагов: " << step_count << endl;

    // Вывод результатов в консоль и запись в CSV
    ofstream out("results.csv");
    if (!out) {
        cerr << "Ошибка открытия файла для записи!" << endl;
        return 1;
    }

    // Заголовок CSV (на английском, чтобы не было проблем с кодировкой в Excel)
    out << "N,t(s),m(kg),P(N),V(m/s),M,Cxa,alpha(deg),theta_c(deg),Cya,omega_z(1/s),theta(deg),y(m),x(m)" << endl;

    // Заголовок для консоли
    cout << "\n==========================================" << endl;
    cout << "Результаты расчетов:" << endl;
    cout << "==========================================" << endl;
    cout << setw(5) << "N"
        << setw(10) << "t, с"
        << setw(10) << "m, кг"
        << setw(12) << "V, м/с"
        << setw(10) << "y, м"
        << setw(10) << "x, м"
        << setw(12) << "alpha, град"
        << endl;
    cout << "------------------------------------------" << endl;

    for (size_t i = 0; i < results.size(); ++i) {
        State& s = results[i];
        AirControl::AtmosphereParams atm = AirControl::AtmosphereCalculator::Calculate(s.y_c);

        double M = s.V / atm.speed_of_sound;
        double Cxa = get_Cxa(M);
        double Cya_alpha = get_Cya(M);
        double alpha_rad = s.theta - s.theta_c;
        double alpha_deg = alpha_rad * 180.0 / M_PI;
        double Cya = Cya_alpha * alpha_rad;
        double P0 = m_dot * W;
        double P = P0 + S_a * (p0N - atm.pressure);

        // Запись в CSV
        out << fixed << setprecision(6)
            << i + 1 << ","
            << s.t << ","
            << s.m << ","
            << P << ","
            << s.V << ","
            << M << ","
            << Cxa << ","
            << alpha_deg << ","
            << s.theta_c * 180.0 / M_PI << ","
            << Cya << ","
            << s.omega_z << ","
            << s.theta * 180.0 / M_PI << ","
            << s.y_c << ","
            << s.x_c << endl;

        // Вывод в консоль каждые 5 шагов
        if (i % 5 == 0 || i == results.size() - 1) {
            cout << fixed << setprecision(2)
                << setw(5) << i + 1
                << setw(10) << s.t
                << setw(10) << s.m
                << setw(12) << s.V
                << setw(10) << s.y_c
                << setw(10) << s.x_c
                << setw(12) << alpha_deg
                << endl;
        }
    }

    out.close();

    // Финальный вывод
    cout << "==========================================" << endl;
    if (!results.empty()) {
        State & final = results.back();
        cout << "\nФинальные параметры:" << endl;
        cout << "Время: " << final.t << " с" << endl;
        cout << "Скорость: " << final.V << " м/с" << endl;
        cout << "Высота: " << final.y_c << " м" << endl;
        cout << "Дальность: " << final.x_c << " м" << endl;
        cout << "Масса: " << final.m << " кг" << endl;
    }

    cout << "\nРезультаты сохранены в файл results.csv" << endl;

    system("pause");

    return 0;
}