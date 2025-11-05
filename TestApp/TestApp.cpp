#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include "Air_control.h"

using namespace std;

// Функция для определения слоя атмосферы
string getAtmosphereLayer(double H) {
    if (H < -2000) return "Вне диапазона ГОСТ";
    else if (H <= 0) return "Нижняя атмосфера";
    else if (H <= 11000) return "Тропосфера";
    else if (H <= 20000) return "Тропопауза";
    else if (H <= 32000) return "Стратосфера (нижняя)";
    else if (H <= 47000) return "Стратосфера (верхняя)";
    else if (H <= 51000) return "Стратопауза";
    else if (H <= 71000) return "Мезосфера (нижняя)";
    else if (H <= 85000) return "Мезосфера (верхняя)";
    else return "Термосфера";
}

// Функция для форматирования чисел
string formatNumber(double value, int precision = 2) {
    stringstream ss;
    ss << fixed << setprecision(precision) << value;
    return ss.str();
}

// Вывод разделителя таблицы
void printSeparator() {
    cout << "+--------------+--------------+------------------+------------------+----------------+--------------+---------------------+" << endl;
}

// Вывод заголовка таблицы
void printTableHeader() {
    printSeparator();
    cout << "| " << left << setw(12) << "Высота (м)" 
         << "| " << setw(12) << "Температура" 
         << "| " << setw(16) << "Давление (Па)" 
         << "| " << setw(16) << "Плотность" 
         << "| " << setw(14) << "Скор. звука" 
         << "| " << setw(12) << "Ускор. g" 
         << "| " << setw(19) << "Слой атмосферы" 
         << "|" << endl;
    printSeparator();
}

// Вывод строки таблицы
void printTableRow(const AirControl::AtmosphereParams& params) {
    string layer = getAtmosphereLayer(params.geopotential_altitude);
    
    cout << "| " << right << setw(12) << formatNumber(params.geometric_altitude, 0)
         << "| " << setw(12) << formatNumber(params.temperature, 2)
         << "| " << setw(16) << formatNumber(params.pressure, 0)
         << "| " << setw(16) << formatNumber(params.density, 4)
         << "| " << setw(14) << formatNumber(params.speed_of_sound, 2)
         << "| " << setw(12) << formatNumber(params.gravity, 5)
         << "| " << left << setw(19) << layer 
         << "|" << endl;
}

// Получение равномерно распределенных высот
vector<double> getUniformAltitudes(double min_alt, double max_alt, int count) {
    vector<double> altitudes;
    if (count <= 1) {
        altitudes.push_back(min_alt);
        return altitudes;
    }
    
    double step = (max_alt - min_alt) / (count - 1);
    
    for (int i = 0; i < count; ++i) {
        altitudes.push_back(min_alt + i * step);
    }
    
    return altitudes;
}

int main() {
    setlocale(LC_ALL, "Russian");
    
    double min_alt, max_alt;
    const int NUM_POINTS = 10;
    
    cout << "=== Air Control ===" << endl;
    cout << "=== ТАБЛИЧНЫЙ ВЫВОД (" << NUM_POINTS << " точек) ===" << endl;
    
    // Ввод данных с проверкой
    while (true) {
        cout << "Введите минимальную высоту (м): ";
        cin >> min_alt;
        cout << "Введите максимальную высоту (м): ";
        cin >> max_alt;
        
        if (min_alt >= max_alt) {
            cout << "Ошибка: минимальная высота должна быть меньше максимальной!" << endl;
            cout << "Пожалуйста, введите значения заново." << endl;
            continue;
        }
        
        if (min_alt < -2000) {
            cout << "Предупреждение: минимальная высота ниже -2000 м может дать неточные результаты!" << endl;
        }
        
        break;
    }
    
    // Получаем равномерно распределенные высоты
    auto altitudes = getUniformAltitudes(min_alt, max_alt, NUM_POINTS);
    
    // Рассчитываем параметры для каждой высоты
    vector<AirControl::AtmosphereParams> results;
    for (double alt : altitudes) {
        results.push_back(AirControl::AtmosphereCalculator::Calculate(alt));
    }
    
    // Выводим основную таблицу
    cout << "\n=== ОСНОВНЫЕ ПАРАМЕТРЫ АТМОСФЕРЫ ===" << endl;
    printTableHeader();
    for (const auto& params : results) {
        printTableRow(params);
    }
    printSeparator();
    
    // Статистика по слоям атмосферы
    cout << "\n=== СТАТИСТИКА ПО СЛОЯМ АТМОСФЕРЫ ===" << endl;
    map<string, int> layerStats;
    for (const auto& params : results) {
        string layer = getAtmosphereLayer(params.geopotential_altitude);
        layerStats[layer]++;
    }
    
    for (const auto& stat : layerStats) {
        cout << stat.first << ": " << stat.second << " точек" << endl;
    }
    
    // Общая информация
    cout << "\n=== ИНФОРМАЦИЯ ===" << endl;
    cout << "Диапазон высот: " << min_alt << " м - " << max_alt << " м" << endl;
    cout << "Количество точек расчета: " << NUM_POINTS << endl;
    cout << "Шаг между точками: " << (max_alt - min_alt) / (NUM_POINTS - 1) << " м" << endl;
    
    // Пример вычисления на уровне моря для проверки
    cout << "\n=== СПРАВОЧНЫЕ ЗНАЧЕНИЯ НА УРОВНЕ МОРЯ ===" << endl;
    auto seaLevel = AirControl::AtmosphereCalculator::Calculate(0);
    cout << "Геопотенциальная высота: " << seaLevel.geopotential_altitude << " м" << endl;
    cout << "Температура: " << seaLevel.temperature << " K" << endl;
    cout << "Давление: " << seaLevel.pressure << " Па" << endl;
    cout << "Плотность: " << seaLevel.density << " кг/м³" << endl;
    cout << "Скорость звука: " << seaLevel.speed_of_sound << " м/с" << endl;
    cout << "Ускорение свободного падения: " << seaLevel.gravity << " м/с²" << endl;
    
    return 0;
}