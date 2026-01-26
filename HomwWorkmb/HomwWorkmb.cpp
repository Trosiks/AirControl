#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <windows.h>
#include <algorithm>
#include "../Air_control/Air_control.h"

using namespace std;

// Определяем константу PI
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Функция для настройки кодировки консоли Windows
void SetConsoleEncoding() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale(LC_ALL, "Russian");
}

// Константы и параметры ЛА
const double m0 = 1255.0;        // начальная масса, кг
const double Jz = 127.0;         // момент инерции, кг*м^2
const double l_d = 0.215;        // запас устойчивости, м
const double S_a = 0.14;         // площадь выходного сечения сопла, м^2
const double S_m = 0.231;        // характерная площадь, м^2
const double g0 = 9.80665;       // ускорение свободного падения на уровне моря, м/с^2
const double p0N = 101325.0;     // нормальное давление у поверхности Земли, Па

// Исходные данные
const double V0 = 65.0;          // начальная скорость, м/с
const double theta_c0_deg = 48.0; // начальный угол наклона траектории, градусы
const double theta_c0 = theta_c0_deg * M_PI / 180.0; // в радианах
const double m_dot = 85.0;       // массовый расход, кг/с
const double W = 2180.0;         // скорость истечения газов, м/с
const double y0_val = 3869.0;    // начальная высота, м
const double omega_z0 = 0.03;    // начальная угловая скорость, с^-1
const double theta0 = theta_c0;  // начальный угол тангажа, рад
const double t0 = 0.0;           // начальное время, с
const double x0 = 0.0;           // начальная координата x, м
const double t_k = 3.72;         // время активного участка, с

// Табличные данные для Cxa(M) и Cya(M)
vector<double> M_values = { 0.01, 0.55, 0.8, 0.9, 1.0, 1.06, 1.1, 1.2, 1.3, 1.4, 2.0, 2.6, 3.4, 6.0, 10.2 };
vector<double> Cxa_values = { 0.30, 0.30, 0.55, 0.70, 0.84, 0.86, 0.87, 0.83, 0.80, 0.79, 0.65, 0.55, 0.50, 0.45, 0.41 };
vector<double> Cya_values = { 0.25, 0.25, 0.25, 0.20, 0.30, 0.31, 0.25, 0.25, 0.25, 0.25, 0.20, 0.25, 0.25, 0.25, 0.15 };

// Структура для хранения состояния системы
struct State {
    double t;       // время, с
    double V;       // скорость, м/с
    double theta_c; // угол наклона траектории, рад
    double x;       // координата x, м
    double y;       // координата y, м
    double omega_z; // угловая скорость, рад/с
    double theta;   // угол тангажа, рад
    double m;       // масса, кг
};

// Структура для хранения результатов записи
struct RecordData {
    double t;
    double m;
    double P;
    double V;
    double M;
    double Cxa;
    double alpha_deg;
    double theta_c_deg;
    double Cya;
    double omega_z;
    double theta_deg;
    double y;
    double x;
};

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

// Функция для вычисления производных (правые части уравнений)
void derivatives(const State& s, State& dsdt, const AirControl::AtmosphereParams& atm, double alpha) {
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
    dsdt.x = s.V * cos(s.theta_c);
    dsdt.y = s.V * sin(s.theta_c);
    dsdt.omega_z = M_z / Jz;
    dsdt.theta = s.omega_z;
    dsdt.t = 1.0;
    dsdt.m = -m_dot;
}

// Метод Эйлера
vector<State> eulerMethod(double dt, double alpha_type, double t_end = 3.8) {
    State s;
    s.t = t0;
    s.V = V0;
    s.theta_c = theta_c0;
    s.x = x0;
    s.y = y0_val;
    s.omega_z = omega_z0;
    s.theta = theta0;
    s.m = m0;

    vector<State> results;
    results.push_back(s);

    for (double t = t0 + dt; t <= t_end + dt / 2; t += dt) {
        // Расчет параметров атмосферы
        AirControl::AtmosphereParams atm = AirControl::AtmosphereCalculator::Calculate(s.y);

        // Вычисление угла атаки
        double alpha;
        if (alpha_type == 0) {
            alpha = 0.0;
        }
        else {
            alpha = s.theta - s.theta_c;
        }

        // Вычисление производных
        State dsdt;
        derivatives(s, dsdt, atm, alpha);

        // Шаг Эйлера
        s.V += dt * dsdt.V;
        s.theta_c += dt * dsdt.theta_c;
        s.x += dt * dsdt.x;
        s.y += dt * dsdt.y;
        s.omega_z += dt * dsdt.omega_z;
        s.theta += dt * dsdt.theta;
        s.t += dt;
        s.m += dt * dsdt.m;

        // Защита от отрицательной массы
        if (s.m < 0) s.m = 0;

        results.push_back(s);
    }

    return results;
}

// Модифицированный метод Эйлера
vector<State> modifiedEulerMethod(double dt, double alpha_type, double t_end = 3.8) {
    State s;
    s.t = t0;
    s.V = V0;
    s.theta_c = theta_c0;
    s.x = x0;
    s.y = y0_val;
    s.omega_z = omega_z0;
    s.theta = theta0;
    s.m = m0;

    vector<State> results;
    results.push_back(s);

    for (double t = t0 + dt; t <= t_end + dt / 2; t += dt) {
        // Расчет параметров атмосферы
        AirControl::AtmosphereParams atm = AirControl::AtmosphereCalculator::Calculate(s.y);

        // Вычисление угла атаки
        double alpha;
        if (alpha_type == 0) {
            alpha = 0.0;
        }
        else {
            alpha = s.theta - s.theta_c;
        }

        // Первая производная
        State k1;
        derivatives(s, k1, atm, alpha);

        // Промежуточное состояние
        State s_temp = s;
        s_temp.V += dt / 2 * k1.V;
        s_temp.theta_c += dt / 2 * k1.theta_c;
        s_temp.x += dt / 2 * k1.x;
        s_temp.y += dt / 2 * k1.y;
        s_temp.omega_z += dt / 2 * k1.omega_z;
        s_temp.theta += dt / 2 * k1.theta;
        s_temp.m += dt / 2 * k1.m;

        // Вторая производная в промежуточной точке
        AirControl::AtmosphereParams atm_temp = AirControl::AtmosphereCalculator::Calculate(s_temp.y);
        State k2;
        derivatives(s_temp, k2, atm_temp, alpha);

        // Полный шаг
        s.V += dt * k2.V;
        s.theta_c += dt * k2.theta_c;
        s.x += dt * k2.x;
        s.y += dt * k2.y;
        s.omega_z += dt * k2.omega_z;
        s.theta += dt * k2.theta;
        s.t += dt;
        s.m += dt * k2.m;

        // Защита от отрицательной массы
        if (s.m < 0) s.m = 0;

        results.push_back(s);
    }

    return results;
}

// Метод Рунге-Кутта 4-го порядка
vector<State> rungeKuttaMethod(double dt, double alpha_type, double t_end = 3.8) {
    State s;
    s.t = t0;
    s.V = V0;
    s.theta_c = theta_c0;
    s.x = x0;
    s.y = y0_val;
    s.omega_z = omega_z0;
    s.theta = theta0;
    s.m = m0;

    vector<State> results;
    results.push_back(s);

    for (double t = t0; t < t_end - dt / 2; t += dt) {
        AirControl::AtmosphereParams atm1 = AirControl::AtmosphereCalculator::Calculate(s.y);
        double alpha1 = (alpha_type == 0) ? 0.0 : (s.theta - s.theta_c);

        State k1;
        derivatives(s, k1, atm1, alpha1);

        State s2 = s;
        s2.V += dt / 2 * k1.V;
        s2.theta_c += dt / 2 * k1.theta_c;
        s2.x += dt / 2 * k1.x;
        s2.y += dt / 2 * k1.y;
        s2.omega_z += dt / 2 * k1.omega_z;
        s2.theta += dt / 2 * k1.theta;
        s2.t += dt / 2;
        s2.m = m0 - m_dot * s2.t;  

        AirControl::AtmosphereParams atm2 = AirControl::AtmosphereCalculator::Calculate(s2.y);
        double alpha2 = (alpha_type == 0) ? 0.0 : (s2.theta - s2.theta_c);

        State k2;
        derivatives(s2, k2, atm2, alpha2);

        State s3 = s;
        s3.V += dt / 2 * k2.V;
        s3.theta_c += dt / 2 * k2.theta_c;
        s3.x += dt / 2 * k2.x;
        s3.y += dt / 2 * k2.y;
        s3.omega_z += dt / 2 * k2.omega_z;
        s3.theta += dt / 2 * k2.theta;
        s3.t += dt / 2;
        s3.m = m0 - m_dot * s3.t;

        AirControl::AtmosphereParams atm3 = AirControl::AtmosphereCalculator::Calculate(s3.y);
        double alpha3 = (alpha_type == 0) ? 0.0 : (s3.theta - s3.theta_c);

        State k3;
        derivatives(s3, k3, atm3, alpha3);

        State s4 = s;
        s4.V += dt * k3.V;
        s4.theta_c += dt * k3.theta_c;
        s4.x += dt * k3.x;
        s4.y += dt * k3.y;
        s4.omega_z += dt * k3.omega_z;
        s4.theta += dt * k3.theta;
        s4.t += dt;
        s4.m = m0 - m_dot * s4.t;

        AirControl::AtmosphereParams atm4 = AirControl::AtmosphereCalculator::Calculate(s4.y);
        double alpha4 = (alpha_type == 0) ? 0.0 : (s4.theta - s4.theta_c);

        State k4;
        derivatives(s4, k4, atm4, alpha4);

        s.V += dt / 6.0 * (k1.V + 2 * k2.V + 2 * k3.V + k4.V);
        s.theta_c += dt / 6.0 * (k1.theta_c + 2 * k2.theta_c + 2 * k3.theta_c + k4.theta_c);
        s.x += dt / 6.0 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
        s.y += dt / 6.0 * (k1.y + 2 * k2.y + 2 * k3.y + k4.y);
        s.omega_z += dt / 6.0 * (k1.omega_z + 2 * k2.omega_z + 2 * k3.omega_z + k4.omega_z);
        s.theta += dt / 6.0 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
        s.t += dt;
        s.m = m0 - m_dot * s.t;  

        if (s.m < 0) s.m = 0;
        results.push_back(s);
    }

    return results;
}
// Функция для интерполяции состояния в заданное время
State interpolateState(const vector<State>& states, double target_time) {
    if (states.empty()) return State();

    // Если время меньше первого, возвращаем первое состояние
    if (target_time <= states[0].t) return states[0];

    // Если время больше последнего, возвращаем последнее состояние
    if (target_time >= states.back().t) return states.back();

    // Ищем два соседних состояния для интерполяции
    for (size_t i = 0; i < states.size() - 1; i++) {
        if (states[i].t <= target_time && target_time <= states[i + 1].t) {
            double t1 = states[i].t;
            double t2 = states[i + 1].t;
            double alpha = (target_time - t1) / (t2 - t1);

            State result;
            result.t = target_time;
            result.V = states[i].V + alpha * (states[i + 1].V - states[i].V);
            result.theta_c = states[i].theta_c + alpha * (states[i + 1].theta_c - states[i].theta_c);
            result.x = states[i].x + alpha * (states[i + 1].x - states[i].x);
            result.y = states[i].y + alpha * (states[i + 1].y - states[i].y);
            result.omega_z = states[i].omega_z + alpha * (states[i + 1].omega_z - states[i].omega_z);
            result.theta = states[i].theta + alpha * (states[i + 1].theta - states[i].theta);
            result.m = states[i].m + alpha * (states[i + 1].m - states[i].m);
            return result;
        }
    }

    return states.back();
}

// Функция для создания записи данных
RecordData createRecord(const State& s, double alpha_type) {
    RecordData record;
    record.t = s.t;
    record.m = s.m;

    // Расчет параметров атмосферы
    AirControl::AtmosphereParams atm = AirControl::AtmosphereCalculator::Calculate(s.y);

    // Тяга
    double P0 = m_dot * W;
    record.P = P0 + S_a * (p0N - atm.pressure);

    record.V = s.V;
    record.M = s.V / atm.speed_of_sound;
    record.Cxa = get_Cxa(record.M);

    // Угол атаки
    double alpha;
    if (alpha_type == 0) {
        alpha = 0.0;
    }
    else {
        alpha = s.theta - s.theta_c;
    }
    record.alpha_deg = alpha * 180.0 / M_PI;

    record.theta_c_deg = s.theta_c * 180.0 / M_PI;

    double Cya_alpha = get_Cya(record.M);
    record.Cya = Cya_alpha * alpha;

    record.omega_z = s.omega_z;
    record.theta_deg = s.theta * 180.0 / M_PI;
    record.y = s.y;
    record.x = s.x;

    return record;
}

// Функция для записи результатов в CSV файл
void writeToCSV(const string& filename, const vector<RecordData>& records) {
    ofstream out(filename);
    if (!out) {
        cerr << "Ошибка открытия файла: " << filename << endl;
        return;
    }

    // Заголовок CSV
    out << "N,t(s),m(kg),P(N),V(m/s),M,Cxa,alpha(deg),theta_c(deg),Cya,omega_z(1/s),theta(deg),y(m),x(m)" << endl;

    for (size_t i = 0; i < records.size(); ++i) {
        const RecordData& r = records[i];
        out << fixed << setprecision(6)
            << i + 1 << ","
            << r.t << ","
            << r.m << ","
            << r.P << ","
            << r.V << ","
            << r.M << ","
            << r.Cxa << ","
            << r.alpha_deg << ","
            << r.theta_c_deg << ","
            << r.Cya << ","
            << r.omega_z << ","
            << r.theta_deg << ","
            << r.y << ","
            << r.x << endl;
    }

    out.close();
    cout << "Результаты записаны в файл: " << filename << endl;
}

// Основная функция
int main() {
    // Настройка кодировки консоли
    SetConsoleEncoding();

    cout << "==========================================" << endl;
    cout << "           Расчет траектории ЛА           " << endl;
    cout << "==========================================" << endl;
    cout << "Исходные данные:" << endl;
    cout << "V0 = " << V0 << " м/с" << endl;
    cout << "theta_c0 = " << theta_c0_deg << " град" << endl;
    cout << "m_dot = " << m_dot << " кг/с" << endl;
    cout << "W = " << W << " м/с" << endl;
    cout << "y0 = " << y0_val << " м" << endl;
    cout << "omega_z0 = " << omega_z0 << " 1/с" << endl;
    cout << "t_k = " << t_k << " с" << endl;
    cout << "==========================================" << endl;

    // Массивы для хранения информации о методах
    vector<string> method_names = { "Euler", "ModifiedEuler", "RungeKutta" };
    vector<double> euler_steps = { 0.1, 0.01, 0.001 };
    vector<double> modified_euler_steps = { 0.1, 0.01 };
    vector<double> runge_kutta_steps = { 0.1 };
    vector<double> alpha_types = { 0, 1 }; // 0: alpha=0, 1: alpha=theta-theta_c

    // Времена для вывода (с шагом 0.1 с и 3.72 с)
    vector<double> output_times;
    for (double t = 0.0; t <= 3.8; t += 0.1) {
        output_times.push_back(t);
    }
    // Добавляем время 3.72 с
    output_times.push_back(3.72);
    // Сортируем времена
    sort(output_times.begin(), output_times.end());

    // Перебор всех методов и углов атаки
    for (size_t method_idx = 0; method_idx < method_names.size(); ++method_idx) {
        string method_name = method_names[method_idx];
        vector<double> steps;

        if (method_name == "Euler") {
            steps = euler_steps;
        }
        else if (method_name == "ModifiedEuler") {
            steps = modified_euler_steps;
        }
        else if (method_name == "RungeKutta") {
            steps = runge_kutta_steps;
        }

        for (double dt : steps) {
            for (double alpha_type : alpha_types) {
                cout << "\nВыполнение расчета: " << method_name
                    << ", шаг = " << dt
                    << ", alpha_type = " << (alpha_type == 0 ? "0" : "theta-theta_c") << endl;

                // Выполнение расчета методом
                vector<State> states;
                if (method_name == "Euler") {
                    states = eulerMethod(dt, alpha_type, 3.8);
                }
                else if (method_name == "ModifiedEuler") {
                    states = modifiedEulerMethod(dt, alpha_type, 3.8);
                }
                else if (method_name == "RungeKutta") {
                    states = rungeKuttaMethod(dt, alpha_type, 3.8);
                }

                // Интерполяция состояний для заданных времен
                vector<RecordData> records;
                for (double t : output_times) {
                    State s_interp = interpolateState(states, t);
                    RecordData record = createRecord(s_interp, alpha_type);
                    records.push_back(record);
                }

                // Создание имени файла
                string filename = method_name + "_step" + to_string(dt).substr(0, 4) +
                    "_alpha" + (alpha_type == 0 ? "0" : "var") + ".csv";

                // Запись в CSV
                writeToCSV(filename, records);

                // Вывод первых нескольких строк для проверки
                cout << "Первые 5 строк результатов:" << endl;
                cout << "t(s)\tm(kg)\tV(m/s)\ty(m)\tx(m)\talpha(deg)" << endl;
                for (int i = 0; i < min(5, (int)records.size()); ++i) {
                    cout << fixed << setprecision(2)
                        << records[i].t << "\t"
                        << records[i].m << "\t"
                        << records[i].V << "\t"
                        << records[i].y << "\t"
                        << records[i].x << "\t"
                        << records[i].alpha_deg << endl;
                }
            }
        }
    }

    cout << "\n==========================================" << endl;
    cout << "Все расчеты завершены!" << endl;
    cout << "Созданы следующие файлы:" << endl;

    // Вывод списка созданных файлов
    vector<string> created_files;
    for (size_t method_idx = 0; method_idx < method_names.size(); ++method_idx) {
        string method_name = method_names[method_idx];
        vector<double> steps;

        if (method_name == "Euler") steps = euler_steps;
        else if (method_name == "ModifiedEuler") steps = modified_euler_steps;
        else if (method_name == "RungeKutta") steps = runge_kutta_steps;

        for (double dt : steps) {
            for (double alpha_type : alpha_types) {
                string filename = method_name + "_step" + to_string(dt).substr(0, 4) +
                    "_alpha" + (alpha_type == 0 ? "0" : "var") + ".csv";
                created_files.push_back(filename);
                cout << "  - " << filename << endl;
            }
        }
    }
    return 0;
}