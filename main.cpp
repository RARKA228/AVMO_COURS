#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <functional>

class Fraction {
    long long num, denom;

    long long computeGCD(long long a, long long b) {
        return b == 0 ? a : computeGCD(b, a % b);
    }

    void reduce() {
        long long divisor = computeGCD(std::abs(num), std::abs(denom));
        num /= divisor;
        denom /= divisor;
        if (denom < 0) {
            num = -num;
            denom = -denom;
        }
    }

public:
    Fraction(long long numerator = 0, long long denominator = 1) : num(numerator), denom(denominator) {
        if (denom == 0) throw std::runtime_error("Denominator cannot be zero");
        reduce();
    }

    long long getNum() const { return num; }
    long long getDenom() const { return denom; }

    double getAbsValue() const {
        return std::abs(static_cast<double>(num) / denom);
    }

    Fraction add(const Fraction& other) const {
        return Fraction(num * other.denom + other.num * denom, denom * other.denom);
    }

    Fraction subtract(const Fraction& other) const {
        return Fraction(num * other.denom - other.num * denom, denom * other.denom);
    }

    Fraction multiply(const Fraction& other) const {
        return Fraction(num * other.num, denom * other.denom);
    }

    Fraction divide(const Fraction& other) const {
        if (other.num == 0) throw std::runtime_error("Division by zero is not allowed");
        return Fraction(num * other.denom, denom * other.num);
    }

    Fraction negate() const {
        return Fraction(-num, denom);
    }

    bool isEqual(const Fraction& other) const {
        return num * other.denom == other.num * denom;
    }

    bool isLessThan(const Fraction& other) const {
        return num * other.denom < other.num * denom;
    }

    Fraction operator+(const Fraction& other) const { return add(other); }
    Fraction operator-(const Fraction& other) const { return subtract(other); }
    Fraction operator*(const Fraction& other) const { return multiply(other); }
    Fraction operator/(const Fraction& other) const { return divide(other); }
    Fraction operator-() const { return negate(); }
    Fraction& operator-=(const Fraction& other) {
        *this = *this - other;
        return *this;
    }
    Fraction& operator+=(const Fraction& other) {
        *this = *this + other;
        return *this;
    }
    bool operator==(const Fraction& other) const { return isEqual(other); }
    bool operator<(const Fraction& other) const { return isLessThan(other); }
    bool operator!=(const Fraction& other) const { return !(*this == other); }
    bool operator>(const Fraction& other) const { return other < *this; }
};

std::ostream& operator<<(std::ostream& os, const Fraction& f) {
    if (f.getDenom() == 1) {
        os << f.getNum();
    }
    else {
        os << f.getNum() << "/" << f.getDenom();
    }
    return os;
}

class Matrix {
    std::vector<std::vector<Fraction>> elements;
    int rowCount, columnCount;

    std::vector<Fraction> objective;
    bool isMaximization;

    std::vector<int> basis;
    std::vector<bool> isArtificialVar;


    long long calculateFactorial(int n) {
        if (n <= 1) return 1;
        return n * calculateFactorial(n - 1);
    }

    long long calculateCombinations(int n, int m) {
        return calculateFactorial(n) / (calculateFactorial(m) * calculateFactorial(n - m));
    }

    std::vector<std::vector<int>> generateAllPossibleBasicSets(int variableCount, int rank) {
        long long expectedCount = calculateCombinations(variableCount, rank);
        std::cout << "\nОжидаемое количество базисных решений: " << expectedCount << "\n";
        std::cout << "Формула расчета: C(" << variableCount << "," << rank << ") = "
            << variableCount << "! / (" << rank << "! * (" << variableCount << "-" << rank << ")!)\n\n";

        std::vector<std::vector<int>> result;
        std::vector<int> combination(rank);

        std::function<void(int, int, int)> generateCombinations =
            [&](int pos, int start, int maxVal) {
            if (pos == rank) {
                result.push_back(combination);
                return;
            }

            for (int i = start; i <= maxVal; i++) {
                combination[pos] = i;
                generateCombinations(pos + 1, i + 1, maxVal);
            }
            };

        generateCombinations(0, 0, variableCount - 1);
        return result;
    }


    bool solveLinearSystem(const std::vector<std::vector<Fraction>>& A, const std::vector<Fraction>& b, std::vector<Fraction>& x) {
        int n = A.size();
        if (n != A[0].size() || n != b.size()) {
            return false;
        }

        std::vector<std::vector<Fraction>> augmented(n, std::vector<Fraction>(n + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][n] = b[i];
        }

        for (int i = 0; i < n; i++) {
            int maxRow = i;
            for (int j = i + 1; j < n; j++) {
                if (augmented[j][i].getAbsValue() > augmented[maxRow][i].getAbsValue()) {
                    maxRow = j;
                }
            }

            if (augmented[maxRow][i] == Fraction(0)) {
                return false;
            }

            if (maxRow != i) {
                for (int j = 0; j <= n; j++) {
                    std::swap(augmented[i][j], augmented[maxRow][j]);
                }
            }


            for (int j = i + 1; j < n; j++) {
                Fraction factor = augmented[j][i] / augmented[i][i];
                for (int k = i; k <= n; k++) {
                    augmented[j][k] = augmented[j][k] - factor * augmented[i][k];
                }
            }
        }


        x.resize(n);
        for (int i = n - 1; i >= 0; i--) {
            Fraction sum = Fraction(0);
            for (int j = i + 1; j < n; j++) {
                sum = sum + augmented[i][j] * x[j];
            }
            x[i] = (augmented[i][n] - sum) / augmented[i][i];
        }

        return true;
    }

    bool isValidSolution(const std::vector<Fraction>& solution) {
        for (int i = 0; i < rowCount; i++) {
            Fraction sum(0);
            for (int j = 0; j < columnCount - 1; j++) {
                sum = sum + elements[i][j] * solution[j];
            }
            if (!(sum == elements[i][columnCount - 1])) {
                return false;
            }
        }

        for (const auto& val : solution) {
            if (val < Fraction(0)) {
                return false;
            }
        }

        return true;
    }

    bool solveWithBasis(const std::vector<int>& basicSet, std::vector<Fraction>& solution) {
        int n = basicSet.size();
        std::vector<std::vector<Fraction>> A(n, std::vector<Fraction>(n));
        std::vector<Fraction> b(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = elements[i][basicSet[j]];
            }
            b[i] = elements[i][columnCount - 1];
        }

        std::vector<Fraction> basicSolution;
        if (!solveLinearSystem(A, b, basicSolution)) {
            return false;
        }

        solution = std::vector<Fraction>(columnCount - 1, Fraction(0));
        for (int i = 0; i < n; i++) {
            solution[basicSet[i]] = basicSolution[i];
        }

        return isValidSolution(solution);
    }


    void updateBasis() {
        basis.clear();
        for (int i = 0; i < rowCount - 1; ++i) {
            for (int j = 0; j < columnCount - 1; ++j) {
                if (elements[i][j] == Fraction(1)) {
                    bool isBasic = true;
                    for (int k = 0; k < rowCount - 1; ++k) {
                        if (k != i && elements[k][j] != Fraction(0)) {
                            isBasic = false;
                            break;
                        }
                    }
                    if (isBasic) {
                        basis.push_back(j);
                        break;
                    }
                }
            }
        }
    }

public:
    Matrix(const std::vector<std::vector<Fraction>>& data, const std::vector<Fraction>& obj, bool isMax)
        : elements(data), objective(obj), isMaximization(isMax) {
        rowCount = elements.size();
        columnCount = rowCount > 0 ? elements[0].size() : 0;
        isArtificialVar.resize(columnCount - 1, false);
        updateBasis();
    }

    void print(const std::string& message = "") const {
        if (!message.empty()) {
            std::cout << message << "\n";
        }
        std::cout << "Матрица " << rowCount << "x" << columnCount << ":\n";
        for (int i = 0; i < rowCount; i++) {
            std::cout << "Строка " << (i + 1) << ": ";
            for (int j = 0; j < columnCount; j++) {
                std::cout << std::setw(10) << elements[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    void jordanGauss() {
    int col = 0;
    for (int row = 0; row < rowCount && col < columnCount - 1; ++row) {
        std::cout << "Шаг " << (row + 1) << ":\n";
        std::cout << "Ищем максимальный элемент в столбце " << (col + 1) << ":\n";
        int leadingRow = row;
        Fraction maxElement(0);
        for (int i = row; i < rowCount - 1; ++i) {
            if (elements[i][col].getAbsValue() > maxElement.getAbsValue()) {
                maxElement = elements[i][col];
                leadingRow = i;
            }
        }
        if (maxElement == Fraction(0)) {
            std::cout << "Нулевой столбец, сдвигаем нулевые строки вниз...\n";
            std::stable_partition(
                elements.begin() + row,
                elements.begin() + (rowCount - 1),
                [](auto &r){
                    return std::any_of(r.begin(), r.end(),
                                       [](const Fraction& f){ return f != Fraction(0); });
                }
            );
            print("После переноса нулевых строк:");
            --row; ++col;
            continue;
        }
        if (leadingRow != row) {
            std::cout << "Меняем строки " << row+1 << " и " << leadingRow+1 << "\n";
            std::swap(elements[row], elements[leadingRow]);
            print();
        }
        Fraction pivot = elements[row][col];
        std::cout << "Делим строку " << row+1 << " на " << pivot << "\n";
        for (int j = col; j < columnCount; ++j) {
            elements[row][j] = elements[row][j] / pivot;
        }
        print();

        for (int i = 0; i < rowCount; ++i) {
            if (i == row) continue;
            Fraction factor = elements[i][col];
            if (factor != Fraction(0)) {
                std::cout << "Вычитаем из строки " << i+1 << " строку "
                          << row+1 << " * " << factor << "\n";
                for (int j = col; j < columnCount; ++j) {
                    elements[i][j] = elements[i][j] - factor * elements[row][j];
                }
                print();
            }
        }


        for (int i = 0; i < rowCount - 1; ++i) {
            bool allZeros = true;
            bool hasNonZero = false;
            for (int j = 0; j < columnCount - 1; ++j) {
                if (elements[i][j] != Fraction(0)) {
                    allZeros = false;
                    hasNonZero = true;
                    break;
                }
            }
            if (allZeros && elements[i][columnCount - 1] != Fraction(0)) {
                throw std::runtime_error("Система противоречива: нет решения");
            }
            if (!allZeros && elements[i][columnCount - 1] < Fraction(0)) {
                bool hasNegativeCoeff = false;
                for (int j = 0; j < columnCount - 1; ++j) {
                    if (elements[i][j] < Fraction(0)) {
                        hasNegativeCoeff = true;
                        break;
                    }
                }
                if (!hasNegativeCoeff) {
                    throw std::runtime_error("Система противоречива: отрицательный свободный член без отрицательных коэффициентов");
                }
            }
        }
        ++col;
    }
    setupObjective();
}


void moveZeroRowsToEnd() {
    std::vector<std::vector<Fraction>> nonZero, zero;
    for (auto &r : elements) {
        bool allZero = std::all_of(r.begin(), r.end(), [](const Fraction& f){ return f == Fraction(0); });
        (allZero ? zero : nonZero).push_back(r);
    }
    elements = nonZero;
    elements.insert(elements.end(), zero.begin(), zero.end());
    rowCount = static_cast<int>(elements.size());
    std::cout << "Матрица после перемещения нулевых строк в конец:\n";
    print();
}


    void checkSolution() {
        bool hasSolution = true;
        std::vector<int> basicVariables;
        std::vector<bool> isBasicVariable(columnCount - 1, false);
        for (int i = 0; i < rowCount; i++) {
            bool allZeros = true;
            for (int j = 0; j < columnCount - 1; j++) {
                if (elements[i][j] != Fraction(0)) {
                    allZeros = false;
                    break;
                }
            }
            if (allZeros && elements[i][columnCount - 1] != Fraction(0)) {
                std::cout << "Система не имеет решений (противоречивое ограничение)\n";
                return;
            }
        }
        for (int i = 0; i < rowCount; i++) {
            int leadingCol = -1;
            for (int j = 0; j < columnCount - 1; j++) {
                if (!(elements[i][j] == Fraction(0))) {
                    if (elements[i][j] == Fraction(1)) {
                        bool isLeading = true;
                        for (int k = 0; k < rowCount; k++) {
                            if (k != i && !(elements[k][j] == Fraction(0))) {
                                isLeading = false;
                                break;
                            }
                        }
                        if (isLeading) {
                            leadingCol = j;
                            isBasicVariable[j] = true;
                            basicVariables.push_back(j);
                            break;
                        }
                    }
                }
            }
        }

        int rank = basicVariables.size();
        int variableCount = columnCount - 1;
        int freeVariables = variableCount - rank;

        if (freeVariables == 0) {
            std::cout << "Система имеет единственное решение\n";
        }
        else {
            std::cout << "Система имеет бесконечно много решений (" << freeVariables << " свободных переменных)\n";
        }

        long long expectedSolutions = calculateCombinations(variableCount, rank);
        std::cout << "\nТеоретическое количество базисных решений: " << expectedSolutions << "\n";

        std::cout << "Поиск базисных решений:\n";
        std::vector<std::vector<int>> allPossibleBasicSets = generateAllPossibleBasicSets(variableCount, rank);
        int validSolutionsCount = 0;

        for (const auto& potentialBasicSet : allPossibleBasicSets) {
            std::vector<Fraction> solution(variableCount, Fraction(0));

            bool validBasis = false;
            try {
                validBasis = solveWithBasis(potentialBasicSet, solution);
            }
            catch (const std::exception& e) {
                validBasis = false;
            }

            std::cout << "Базис { ";
            for (int idx : potentialBasicSet) {
                std::cout << "x" << idx + 1 << " ";
            }
            std::cout << "} ";

            if (validBasis) {
                validSolutionsCount++;
                std::cout << "является допустимым базисным решением.\n";
                for (int j = 0; j < variableCount; j++) {
                    std::cout << "x" << j + 1 << " = " << solution[j] << "\n";
                }
                std::cout << "\n";
            }
            else {
                std::cout << "не является допустимым базисным решением.\n";
            }
        }

        std::cout << "Найдено базисных решений: " << validSolutionsCount << "\n";
        if (validSolutionsCount != expectedSolutions) {
            std::cout << "Количество найденных базисных решений (" << validSolutionsCount
                << ") отличается от теоретического (" << expectedSolutions << ")\n";
        }
    }


    void prepareForSimplex() {
        for (int i = 0; i < rowCount - 1; ++i) {
            if (elements[i][columnCount - 1] < Fraction(0)) {
                bool canResolve = false;
                for (int j = 0; j < columnCount - 1; ++j) {
                    if (elements[i][j] < Fraction(0)) {
                        canResolve = true;
                        break;
                    }
                }
                if (!canResolve) {
                    throw std::runtime_error("Система противоречива: отрицательный свободный член и нет отрицательных коэффициентов");
                }
                for (int j = 0; j < columnCount; ++j) {
                    elements[i][j] = -elements[i][j];
                }
            }
        }
    }

    void printSimplexTable(int currentPivotColumn = -1, const std::string& message = "") {
        updateBasis();
        if (!message.empty()) {
            std::cout << message << "\n";
        }
        std::cout << "Бп  1        ";
        for (int j = 0; j < columnCount - 1; ++j) {
            std::cout << "x" << j+1 << std::setw(10);
        }
        std::cout << "\n";
        for (int i = 0; i < rowCount - 1; ++i) {
            bool isBasic = false;
            int basicVar = -1;
            for (int j = 0; j < columnCount - 1; ++j) {
                if (elements[i][j] == Fraction(1)) {
                    bool isLeading = true;
                    for (int k = 0; k < rowCount - 1; ++k) {
                        if (k != i && elements[k][j] != Fraction(0)) {
                            isLeading = false;
                            break;
                        }
                    }
                    if (isLeading) {
                        basicVar = j;
                        isBasic = true;
                        break;
                    }
                }
            }
            if (isBasic) {
                std::cout << "x" << basicVar + 1 << "  ";
            } else {
                std::cout << "-   "; // Если нет базисной переменной
            }
            std::cout << elements[i][columnCount - 1] << "";
            for (int j = 0; j < columnCount - 1; ++j) {
                std::cout << std::setw(10) << elements[i][j] << "";
            }
            if (currentPivotColumn >= 0 && elements[i][currentPivotColumn] > Fraction(0)) {
                Fraction ratio = elements[i][columnCount - 1] / elements[i][currentPivotColumn];
                std::cout << "  (Ratio: " << ratio << ")";
            }
            std::cout << "\n";
        }

        // Вывод Z-строки
        std::cout << "Z   " << elements.back()[columnCount - 1] << "";
        for (int j = 0; j < columnCount - 1; ++j) {
            std::cout << std::setw(10) << elements.back()[j] << "";
        }
        std::cout << "\n\n";
    }


    void setupObjective() {
        std::vector<Fraction> zRow(columnCount, Fraction(0));
        for (int j = 0; j < columnCount - 1; ++j) {
            zRow[j] = -objective[j];
        }
        for (int i = 0; i < rowCount; ++i) {
            for (int j = 0; j < columnCount - 1; ++j) {
                if (elements[i][j] == Fraction(1)) {
                    bool isBasic = true;
                    for (int k = 0; k < rowCount; ++k) {
                        if (k != i && elements[k][j] != Fraction(0)) {
                            isBasic = false;
                            break;
                        }
                    }
                    if (isBasic) {
                        Fraction coeff = zRow[j];
                        for (int k = 0; k < columnCount; ++k) {
                            zRow[k] = zRow[k] - coeff * elements[i][k];
                        }
                    }
                }
            }
        }
        if (isMaximization) {
            for (auto & v : zRow) v = -v;
        }
        elements.back() = zRow;
    }


  void simplexMethod() {
    try {
        prepareForSimplex();
    } catch (const std::runtime_error& e) {
        std::cout << "Система не имеет решений: " << e.what() << "\n";
        return;
    }
    setupObjective();
    std::cout << "\nНачальная симплекс-таблица:\n";
    printSimplexTable();
    int iteration = 1;

    while (true) {

        std::cout << "\n=== Итерация " << iteration << " — текущее состояние таблицы ===\n";
        printSimplexTable(-1, "Таблица перед поиском разрешающего столбца:");

        for (int i = 0; i < rowCount - 1; ++i) {
            if (elements[i][columnCount - 1] < Fraction(0)) {
                bool hasNegativeCoeff = false;
                for (int j = 0; j < columnCount - 1; ++j) {
                    if (elements[i][j] < Fraction(0)) {
                        hasNegativeCoeff = true;
                        break;
                    }
                }
                if (!hasNegativeCoeff) {
                    std::cout << "Система противоречива: нет решения\n";
                    printCurrentSolution();
                    return;
                }
            }
        }


        int pivotCol = -1;
        Fraction maxVal(0);
        for (int j = 0; j < columnCount - 1; ++j) {
            Fraction current = elements.back()[j];
            if (isMaximization ? (current < 0) : (current > 0)) {
                if (pivotCol == -1 || current.getAbsValue() > maxVal.getAbsValue()) {
                    pivotCol = j;
                    maxVal = current;
                }
            }
        }

        if (pivotCol == -1) {

            bool isFeasible = true;
            for (int i = 0; i < rowCount - 1; ++i) {
                if (elements[i][columnCount - 1] < Fraction(0)) {
                    isFeasible = false;
                    break;
                }
            }
            if (!isFeasible) {
                std::cout << "Система не имеет допустимых решений: отрицательные значения в столбце свободных членов\n";
                printCurrentSolution();
                return;
            }

            bool hasArtificialInZ = false;
            for (int j = 0; j < columnCount - 1; ++j) {
                if (isArtificialVar[j] && elements.back()[j] != Fraction(0)) {
                    hasArtificialInZ = true;
                    break;
                }
            }

            if (hasArtificialInZ) {
                std::cout << "Задача не имеет допустимых решений!\n";
                return;
            }

            std::cout << "\nВсе коэффициенты в Z-строке "
                      << (isMaximization ? "неотрицательны" : "неположительны")
                      << ". Оптимальное решение найдено.\n";
            break;
        }

        std::cout << "\nШаг " << iteration << ":\n";
        std::cout << "Выбираем разрешающий столбец x" << pivotCol + 1
                  << " (коэффициент: " << elements.back()[pivotCol] << ")\n";


        bool hasPositive = false;
        for (int i = 0; i < rowCount - 1; ++i) {
            if (elements[i][pivotCol] > Fraction(0)) {
                hasPositive = true;
                break;
            }
        }
        if (!hasPositive) {
            std::cout << "Нет положительных элементов в столбце. Функция не ограничена!\n";
            printCurrentSolution();
            return;
        }

        int pivotRow = -1;
        Fraction minRatio;
        bool first = true;
        std::vector<Fraction> ratios(rowCount - 1);

        std::cout << "Симплексные отношения:\n";
        for (int i = 0; i < rowCount - 1; ++i) {
            if (elements[i][pivotCol] > 0) {
                ratios[i] = elements[i][columnCount - 1] / elements[i][pivotCol];
                std::cout << "Строка " << i+1 << ": " << elements[i][columnCount - 1]
                          << " / " << elements[i][pivotCol] << " = " << ratios[i] << "\n";

                if (first || ratios[i] < minRatio) {
                    minRatio = ratios[i];
                    pivotRow = i;
                    first = false;
                }
            } else {
                ratios[i] = Fraction(-1);
                std::cout << "Строка " << i+1 << ": не рассматривается (элемент ≤ 0)\n";
            }
        }

        std::cout << "Выбираем разрешающую строку " << pivotRow + 1
                  << " с минимальным отношением " << minRatio << "\n";

        Fraction pivot = elements[pivotRow][pivotCol];
        std::cout << "\nШаг " << iteration << ".1: Делим строку " << pivotRow + 1
                  << " на ведущий элемент " << pivot << ":\n";

        for (int j = 0; j < columnCount; ++j) {
            elements[pivotRow][j] = elements[pivotRow][j] / pivot;
        }

        printSimplexTable(pivotCol, "Таблица после выбора разрешающего столбца");
        std::cout << "После итерации " << iteration << " получилась такая таблица:\n";
        printSimplexTable();

        std::cout << "Шаг " << iteration << ".2: Обнуляем столбец x" << pivotCol + 1 << ":\n";
        for (int i = 0; i < rowCount; ++i) {
            if (i == pivotRow) continue;

            Fraction factor = elements[i][pivotCol];
            if (factor != 0) {
                std::cout << "Строка " << i+1 << " = Строка " << i+1 << " - ("
                          << factor << " * Строка " << pivotRow + 1 << ")\n";

                for (int j = 0; j < columnCount; ++j) {
                    elements[i][j] = elements[i][j] - factor * elements[pivotRow][j];
                }
            }
        }
        printSimplexTable(-1, "Таблица после преобразований");

        basis[pivotRow] = pivotCol;

        std::cout << "Текущее значение Z = "
                  << elements.back()[columnCount - 1]
                  << (isMaximization ? " (неоптимальное)\n" : "\n");

        basis[pivotRow] = pivotCol;
        iteration++;
    }

    printFinalSolution();
}
    void printCurrentSolution() {
        std::cout << "\nТекущее решение:\n";
        bool isFeasible = true;
        for (int j = 0; j < columnCount - 1; ++j) {
            bool isBasic = false;
            for (size_t i = 0; i < basis.size(); ++i) {
                if (basis[i] == j) {
                    Fraction value = elements[i][columnCount - 1];
                    if (value < Fraction(0)) {
                        isFeasible = false;
                    }
                    std::cout << "x" << j + 1 << " = " << value << "\n";
                    isBasic = true;
                    break;
                }
            }
            if (!isBasic) {
                std::cout << "x" << j + 1 << " = 0\n";
            }
        }
        if (!isFeasible) {
            std::cout << "Решение недопустимо (отрицательные значения переменных)!\n";
        }
        std::cout << "Z = " << elements.back()[columnCount - 1] << "\n";
    }
    void printFinalSolution() {
        std::cout << "\nФинальное решение:\n";
        printCurrentSolution();

        bool hasAlternative = false;
        for (int j = 0; j < columnCount - 1; ++j) {
            if (elements.back()[j] == 0) {
                hasAlternative = true;
                break;
            }
        }

        if (hasAlternative) {
            std::cout << "\n";
        }
    }

};

int main() {
    std::ifstream file("input.txt");
    if (!file) {
        std::cerr << "Ошибка: Невозможно открыть файл input.txt\n";
        return 1;
    }

    int rows, columns;
    file >> rows >> columns;

    std::vector<std::vector<Fraction>> fullCoeffs(rows, std::vector<Fraction>(columns));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            long long num;
            if (!(file >> num)) {
                std::cerr << "Ошибка: Неверные данные во входном файле\n";
                return 1;
            }
            fullCoeffs[i][j] = Fraction(num);
        }
    }


    bool isMax;
    if (!(file >> isMax)) {
        std::cerr << "Ошибка: Не удалось прочитать флаг isMax\n";
        return 1;
    }

    std::vector<Fraction> objectiveCoefficients;
    objectiveCoefficients.reserve(columns - 1);
    for (int j = 0; j < columns - 1; ++j) {
        objectiveCoefficients.push_back(fullCoeffs.back()[j]);
    }

    Matrix matrix(fullCoeffs, objectiveCoefficients, isMax);

    std::cout << "Начальная матрица:\n";
    matrix.print();

    matrix.moveZeroRowsToEnd();
    matrix.jordanGauss();
    matrix.simplexMethod();
    matrix.checkSolution();

    return 0;
}
