#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include "Air_control.h"

using namespace std;

// Функция для форматирования чисел
string formatNumber(double value, int precision = 2) {
    stringstream ss;
    ss << fixed << setprecision(precision) << value;
    return ss.str();
}

// Функция для вывода разделительной линии
void printSeparator(int width) {
    cout << "+" << string(width - 2, '-') << "+" << endl;
}

// Функция для центрирования текста
string centerText(const string& text, int width) {
    if (text.length() >= width) return text;
    int padding = (width - text.length()) / 2;
    return string(padding, ' ') + text + string(width - text.length() - padding, ' ');
}

// Функция для вывода строки таблицы
void printTableRow(const string& label, const string& value, const string& unit, int width) {
    string fullText = label + ": " + value + " " + unit;
    cout << "| " << left << setw(width - 4) << fullText << " |" << endl;
}

// Функция для вывода параметров на конкретной высоте
void printAltitudeParameters(double altitude) {
    const int TABLE_WIDTH = 60;

    // Рассчитываем параметры
    auto params = AirControl::AtmosphereCalculator::Calculate(altitude);

    // Выводим заголовок
    cout << "\n";
    printSeparator(TABLE_WIDTH);
    string title = "ПАРАМЕТРЫ АТМОСФЕРЫ НА ВЫСОТЕ " + formatNumber(altitude, 0) + " М";
    cout << "|" << centerText(title, TABLE_WIDTH - 2) << "|" << endl;
    printSeparator(TABLE_WIDTH);

    // Выводим параметры
    printTableRow("Геопотенциальная высота", formatNumber(params.geopotential_altitude, 0), "м", TABLE_WIDTH);
    printTableRow("Геометрическая высота", formatNumber(params.geometric_altitude, 0), "м", TABLE_WIDTH);
    printTableRow("Температура", formatNumber(params.temperature, 2), "K", TABLE_WIDTH);
    printTableRow("Давление", formatNumber(params.pressure, 0), "Па", TABLE_WIDTH);
    printTableRow("Плотность", formatNumber(params.density, 4), "кг/м³", TABLE_WIDTH);
    printTableRow("Скорость звука", formatNumber(params.speed_of_sound, 2), "м/с", TABLE_WIDTH);
    printTableRow("Ускорение своб. падения", formatNumber(params.gravity, 5), "м/с²", TABLE_WIDTH);

    printSeparator(TABLE_WIDTH);
}

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

// Функция для вывода информации о слое атмосферы
void printLayerInfo(double altitude) {
    auto params = AirControl::AtmosphereCalculator::Calculate(altitude);
    string layer = getAtmosphereLayer(params.geopotential_altitude);

    cout << "\nСлой атмосферы: " << layer << endl;
}

int main() {
    setlocale(LC_ALL, "Russian");

    double altitude;

    cout << "==========================================" << endl;
    cout << "      Параметры на заданной высоте        " << endl;
    cout << "==========================================" << endl;
    cout << "              ГОСТ 4401-81                " << endl;
    cout << "Диапазон высот: от -2000 м до 94000 м" << endl;
    cout << endl;

    // Ввод высоты с проверкой
    while (true) {
        cout << "Введите высоту (м): ";
        cin >> altitude;

        if (altitude < -2000 || altitude > 94000) {
            cout << "Ошибка: высота должна быть в диапазоне от -2000 до 94000 м!" << endl;
            cout << "Пожалуйста, введите корректное значение." << endl;
            continue;
        }

        break;
    }

    // Вывод параметров
    printAltitudeParameters(altitude);
    printLayerInfo(altitude);

    // Дополнительная информация
    cout << "\nДля сравнения - параметры на уровне моря:" << endl;
    auto seaLevel = AirControl::AtmosphereCalculator::Calculate(0);
    cout << "Температура: " << formatNumber(seaLevel.temperature, 2) << " K" << endl;
    cout << "Давление: " << formatNumber(seaLevel.pressure, 0) << " Па" << endl;
    cout << "Плотность: " << formatNumber(seaLevel.density, 4) << " кг/м³" << endl;
    cout << "Скорость звука: " << formatNumber(seaLevel.speed_of_sound, 2) << " м/с" << endl;

    return 0;
}