#include <windows.h>
#include <windowsx.h>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "Air_control.h"

using namespace std;

// Глобальные переменные
HWND g_hMainWindow;
vector<AirControl::AtmosphereParams> g_data;
double g_minAlt = -2000;
double g_maxAlt = 94000;
int g_selectedParameter = 0;

// Константы с увеличенными отступами
const int PANEL_WIDTH = 200;      // Ширина панели с кнопками
const int MARGIN_LEFT = 80;       // Отступ слева для подписей значений
const int MARGIN_TOP = 60;        // Отступ сверху для заголовка
const int MARGIN_RIGHT = 20;      // Отступ справа
const int MARGIN_BOTTOM = 80;     // Отступ снизу для подписей высоты
const int MARGIN_BUTTONS = 20;    // Отступ для кнопок

// Размеры графика (будут вычисляться динамически)
int g_chartWidth = 600;
int g_chartHeight = 400;
int g_chartLeft = MARGIN_LEFT;
int g_chartTop = MARGIN_TOP;

// Функция для расчета данных
void CalculateData() {
    g_data = AirControl::AtmosphereCalculator::CalculateRange(g_minAlt, g_maxAlt, 1000);
}

// Функция для получения значения выбранного параметра
double GetParameterValue(const AirControl::AtmosphereParams& params, int paramIndex) {
    switch (paramIndex) {
    case 0: return params.temperature;
    case 1: return params.pressure;
    case 2: return params.density;
    case 3: return params.speed_of_sound;
    case 4: return params.gravity;
    default: return params.temperature;
    }
}

// Функция для получения имени параметра
wstring GetParameterName(int paramIndex) {
    switch (paramIndex) {
    case 0: return L"Температура";
    case 1: return L"Давление";
    case 2: return L"Плотность";
    case 3: return L"Скорость звука";
    case 4: return L"Ускорение";
    default: return L"Температура";
    }
}

// Функция для получения единицы измерения
wstring GetParameterUnit(int paramIndex) {
    switch (paramIndex) {
    case 0: return L"K";
    case 1: return L"Па";
    case 2: return L"кг/м³";
    case 3: return L"м/с";
    case 4: return L"м/с²";
    default: return L"";
    }
}

// Функция для отрисовки текста
void DrawText(HDC hdc, int x, int y, const wstring& text) {
    TextOutW(hdc, x, y, text.c_str(), (int)text.length());
}

// Функция для вычисления размеров при изменении окна
void CalculateLayout(RECT& clientRect) {
    g_chartWidth = clientRect.right - clientRect.left - PANEL_WIDTH - MARGIN_LEFT - MARGIN_RIGHT;
    g_chartHeight = clientRect.bottom - clientRect.top - MARGIN_TOP - MARGIN_BOTTOM;
    g_chartLeft = MARGIN_LEFT;
    g_chartTop = MARGIN_TOP;
}

// Функция для отрисовки сетки и осей
void DrawGrid(HDC hdc, RECT& clientRect) {
    HPEN gridPen = CreatePen(PS_SOLID, 1, RGB(200, 200, 200));
    HPEN axisPen = CreatePen(PS_SOLID, 2, RGB(0, 0, 0));
    HPEN oldPen = (HPEN)SelectObject(hdc, gridPen);

    // Горизонтальные линии сетки (высота)
    for (int i = 0; i <= 10; i++) {
        int y = g_chartTop + i * g_chartHeight / 10;
        MoveToEx(hdc, g_chartLeft, y, NULL);
        LineTo(hdc, g_chartLeft + g_chartWidth, y);

        // Подписи высоты 
        double alt = g_maxAlt - i * (g_maxAlt - g_minAlt) / 10; 
        wstring label = to_wstring((int)alt) + L" м";
        DrawText(hdc, g_chartLeft - 70, y - 8, label);
    }

    // Вертикальные линии сетки 
    for (int i = 0; i <= 10; i++) {
        int x = g_chartLeft + i * g_chartWidth / 10;
        MoveToEx(hdc, x, g_chartTop, NULL);
        LineTo(hdc, x, g_chartTop + g_chartHeight);
    }

    // Оси
    SelectObject(hdc, axisPen);

    // Ось Y 
    MoveToEx(hdc, g_chartLeft, g_chartTop, NULL);
    LineTo(hdc, g_chartLeft, g_chartTop + g_chartHeight);

    // Ось X 
    MoveToEx(hdc, g_chartLeft, g_chartTop + g_chartHeight, NULL);
    LineTo(hdc, g_chartLeft + g_chartWidth, g_chartTop + g_chartHeight);

    SelectObject(hdc, oldPen);
    DeleteObject(gridPen);
    DeleteObject(axisPen);
}

// Функция для отрисовки графика
void DrawChart(HDC hdc, RECT& clientRect) {
    if (g_data.empty()) return;

    // Находим минимальное и максимальное значение выбранного параметра
    double minVal = GetParameterValue(g_data[0], g_selectedParameter);
    double maxVal = minVal;

    for (const auto& point : g_data) {
        double val = GetParameterValue(point, g_selectedParameter);
        if (val < minVal) minVal = val;
        if (val > maxVal) maxVal = val;
    }

    // Добавляем немного места
    double range = maxVal - minVal;
    if (range > 0) {
        minVal -= range * 0.05;
        maxVal += range * 0.05;
    }
    else {
        minVal -= 1;
        maxVal += 1;
    }

    // Создаем перо для графика
    HPEN chartPen = CreatePen(PS_SOLID, 3, RGB(0, 100, 200));
    HPEN oldPen = (HPEN)SelectObject(hdc, chartPen);

    // Рисуем график (теперь X - значение параметра, Y - высота)
    for (size_t i = 0; i < g_data.size() - 1; i++) {
        // X координата - значение параметра
        double x1 = g_chartLeft + (GetParameterValue(g_data[i], g_selectedParameter) - minVal) * g_chartWidth / (maxVal - minVal);
        // Y координата - высота (инвертируем для правильного направления)
        double y1 = g_chartTop + g_chartHeight - (g_data[i].geopotential_altitude - g_minAlt) * g_chartHeight / (g_maxAlt - g_minAlt);

        double x2 = g_chartLeft + (GetParameterValue(g_data[i + 1], g_selectedParameter) - minVal) * g_chartWidth / (maxVal - minVal);
        double y2 = g_chartTop + g_chartHeight - (g_data[i + 1].geopotential_altitude - g_minAlt) * g_chartHeight / (g_maxAlt - g_minAlt);

        MoveToEx(hdc, (int)x1, (int)y1, NULL);
        LineTo(hdc, (int)x2, (int)y2);
    }

    // Подписи осей
    wstring xLabel = GetParameterName(g_selectedParameter) + L", " + GetParameterUnit(g_selectedParameter);
    DrawText(hdc, g_chartLeft + g_chartWidth / 2 - 80, g_chartTop + g_chartHeight + 40, xLabel);

    // Подпись оси Y 
    DrawText(hdc, 10, g_chartTop + g_chartHeight / 2 - 50, L"Высота, м");

    // Значения на оси X
    for (int i = 0; i <= 5; i++) {
        double val = minVal + i * (maxVal - minVal) / 5;
        wstringstream ss;
        ss << fixed << setprecision(2) << val;
        wstring label = ss.str();

        int x = g_chartLeft + i * g_chartWidth / 5;
        DrawText(hdc, x - 15, g_chartTop + g_chartHeight + 20, label);
    }

    SelectObject(hdc, oldPen);
    DeleteObject(chartPen);
}

// Функция для обновления позиций кнопок
void UpdateButtonPositions(RECT& clientRect) {
    int buttonWidth = PANEL_WIDTH - MARGIN_BUTTONS * 2;
    int buttonHeight = 30;
    int buttonLeft = clientRect.right - PANEL_WIDTH + MARGIN_BUTTONS;
    int buttonTop = MARGIN_TOP;
    int buttonSpacing = 10;

    // Заголовок панели
    HWND hTitle = GetDlgItem(g_hMainWindow, 200);
    if (!hTitle) {
        CreateWindowA("STATIC", "Параметры:", WS_VISIBLE | WS_CHILD | SS_LEFT,
            buttonLeft, buttonTop, buttonWidth, 20, g_hMainWindow, (HMENU)200, NULL, NULL);
    }
    else {
        SetWindowPos(hTitle, NULL, buttonLeft, buttonTop, buttonWidth, 20, SWP_NOZORDER);
    }

    // Кнопки параметров
    for (int i = 100; i <= 104; i++) {
        HWND hButton = GetDlgItem(g_hMainWindow, i);
        if (hButton) {
            SetWindowPos(hButton, NULL,
                buttonLeft,
                buttonTop + 30 + (i - 100) * (buttonHeight + buttonSpacing),
                buttonWidth,
                buttonHeight,
                SWP_NOZORDER);
        }
    }
}

// Функция обработки сообщений окна
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    switch (uMsg) {
    case WM_CREATE: {
        // Создаем элементы управления
        CreateWindowA("BUTTON", "Температура", WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
            0, 0, 100, 30, hwnd, (HMENU)100, NULL, NULL);
        CreateWindowA("BUTTON", "Давление", WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
            0, 0, 100, 30, hwnd, (HMENU)101, NULL, NULL);
        CreateWindowA("BUTTON", "Плотность", WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
            0, 0, 100, 30, hwnd, (HMENU)102, NULL, NULL);
        CreateWindowA("BUTTON", "Скорость звука", WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
            0, 0, 100, 30, hwnd, (HMENU)103, NULL, NULL);
        CreateWindowA("BUTTON", "Ускорение", WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
            0, 0, 100, 30, hwnd, (HMENU)104, NULL, NULL);

        // Устанавливаем первую кнопку как выбранную
        CheckRadioButton(hwnd, 100, 104, 100);

        // Рассчитываем данные при создании окна
        CalculateData();
        return 0;
    }

    case WM_SIZE: {
        RECT clientRect;
        GetClientRect(hwnd, &clientRect);
        CalculateLayout(clientRect);
        UpdateButtonPositions(clientRect);
        InvalidateRect(hwnd, NULL, TRUE); // Перерисовываем окно
        return 0;
    }

    case WM_COMMAND: {
        int id = LOWORD(wParam);
        if (id >= 100 && id <= 104) {
            g_selectedParameter = id - 100;
            CheckRadioButton(hwnd, 100, 104, id);
            InvalidateRect(hwnd, NULL, TRUE); // Перерисовываем окно
        }
        return 0;
    }

    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);

        RECT clientRect;
        GetClientRect(hwnd, &clientRect);

        // Очищаем фон
        HBRUSH background = CreateSolidBrush(RGB(255, 255, 255));
        FillRect(hdc, &clientRect, background);
        DeleteObject(background);

        // Рисуем разделительную линию между графиком и панелью
        HPEN separatorPen = CreatePen(PS_SOLID, 2, RGB(200, 200, 200));
        HPEN oldPen = (HPEN)SelectObject(hdc, separatorPen);
        int separatorX = clientRect.right - PANEL_WIDTH;
        MoveToEx(hdc, separatorX, clientRect.top, NULL);
        LineTo(hdc, separatorX, clientRect.bottom);
        SelectObject(hdc, oldPen);
        DeleteObject(separatorPen);

        // Рисуем сетку и график
        DrawGrid(hdc, clientRect);
        DrawChart(hdc, clientRect);

        EndPaint(hwnd, &ps);
        return 0;
    }

    case WM_DESTROY: {
        PostQuitMessage(0);
        return 0;
    }

    default:
        return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }
}

// Точка входа приложения
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    // Регистрируем класс окна
    WNDCLASS wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = L"AtmosphereGraphApp";
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);

    RegisterClass(&wc);

    // Создаем окно (увеличиваем начальный размер для лучшего отображения)
    g_hMainWindow = CreateWindow(
        L"AtmosphereGraphApp",
        L"Air Control - Визуализация параметров атмосферы по ГОСТ 4401-81",
        WS_OVERLAPPEDWINDOW | WS_SIZEBOX,
        CW_USEDEFAULT, CW_USEDEFAULT, 1200, 800,  // Увеличили начальный размер
        NULL, NULL, hInstance, NULL
    );

    if (g_hMainWindow == NULL) {
        return 0;
    }

    // Показываем окно
    ShowWindow(g_hMainWindow, nCmdShow);
    UpdateWindow(g_hMainWindow);

    // Цикл сообщений
    MSG msg = {};
    while (GetMessage(&msg, NULL, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return 0;
}