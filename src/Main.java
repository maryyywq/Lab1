import java.io.*;
import java.util.*;

public class Main {
    private static final double EPSILON = 1e-5; //Константа для точности вычислений
    private static final int PRECISION = (int) Math.abs(Math.log10(EPSILON));; //Количество знаков после запятой в выводе

    //Класс для хранения результата работы
    public static class Result {
        public double[] solution;   //Найденное решение
        public double[][] inverseMatrix; //Обратная матрица
        public double determinant;      //Определитель
        public double[] residuals; //Вектор невязок

        public Result(double[] solution, double[][] inverseMatrix, double determinant, double[] residuals) {
            this.solution = solution;
            this.determinant = determinant;
            this.inverseMatrix = inverseMatrix;
            this.residuals = residuals;
        }
    }

    public static void main(String[] args) {
        try {
            double[][] A; //Коэффициенты системы уравнений
            double[] b; //Вектор правой части

            //Считывание матрицы из файла
            File inputFile = new File("input.txt");
            Scanner scanner = new Scanner(inputFile);
            int n = scanner.nextInt(); //Читаем размерность матрицы
            A = new double[n][n];
            b = new double[n];

            //Заполняем матрицу A и вектор b
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    A[i][j] = scanner.nextDouble();
                }
                b[i] = scanner.nextDouble();
            }
            scanner.close();

            //Решение системы методом Гаусса с выбором главного элемента
            Result res = gaussianElimination(A, b);

            //Запись и вывод результата
            BufferedWriter writer = new BufferedWriter(new FileWriter("output.txt"));
            if (Double.isNaN(res.determinant)) {
                writer.write("Матрица системы вырожденная. Решение невозможно.\n");
                System.out.println("Матрица системы вырожденная. Решение невозможно.");
            } else {
                writer.write("Решение системы:\n");
                System.out.println("Решение системы:");
                for (double x : res.solution) {
                    writer.write(String.format("%." + PRECISION + "f\n", x));
                    System.out.println(String.format("%." + PRECISION + "f", x));
                }
                writer.write("\nНевязки:\n");
                System.out.println("\nНевязки:");
                for (double r : res.residuals) {
                    writer.write(String.format("%." + PRECISION + "f\n", r));
                    System.out.println(String.format("%." + PRECISION + "f", r));
                }
                writer.write("\nОпределитель: " + String.format("%." + PRECISION + "f\n", res.determinant));
                System.out.println("\nОпределитель: " + String.format("%." + PRECISION + "f", res.determinant));
                writer.write("\nОбратная матрица:\n");
                System.out.println("\nОбратная матрица:");
                for (double[] row : res.inverseMatrix) {
                    for (double value : row) {
                        writer.write(String.format("%." + PRECISION + "f ", value));
                        System.out.print(String.format("%." + PRECISION + "f ", value));
                    }
                    writer.write("\n");
                    System.out.println();
                }
            }
            writer.close();
        } catch (IOException e) {
            System.out.println("Ошибка при работе с файлами: " + e.getMessage());
        }
    }

    public static Result gaussianElimination(double[][] A, double[] b) {
        int n = A.length;
        double[][] identityMatrix = new double[n][n]; //Единичная матрица (E)
        for (int i = 0; i < n; i++) {
            identityMatrix[i][i] = 1; //Ставим 1 по главной диагонали в единичной матрице
        }
        int swapCount = 0; //Количество перестановок строк

        double det = 1.0; //Определитель
        double[] solution = new double[n];
        double[][] inverseMatrix = new double[n][n];
        double[] residuals = new double[n];

        double[][] oldA = new double[n][n]; //Матрица для хранения старых коэффициентов, чтобы посчитать невязки
        for (int i = 0; i < n; i++) System.arraycopy(A[i], 0, oldA[i], 0, n); //Скопируем все элементы
        double[] oldB = new double[n]; //Вектор свободных членов со старыми коэффициентами, чтобы посчитать невязки
        System.arraycopy(b, 0, oldB, 0, n); //Скопируем все элементы


        //Прямой ход метода Гаусса с выбором главного элемента
        for (int i = 0; i < n; i++) {
            //Выбор строки с наибольшим (главным) элементом
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }

            //Проверка на вырожденность матрицы
            if (Math.abs(A[maxRow][i]) < EPSILON) {
                return new Result(null, null, Double.NaN, null);
            }

            //Перестановка строк, если необходимо
            if (i != maxRow) {
                //Перестановка в главной матрице
                double[] temp = A[i];
                A[i] = A[maxRow];
                A[maxRow] = temp;

                //В единичной
                temp = identityMatrix[i];
                identityMatrix[i] = identityMatrix[maxRow];
                identityMatrix[maxRow] = temp;

                //В столбце свободных членов
                double t = b[i];
                b[i] = b[maxRow];
                b[maxRow] = t;

                swapCount++;
            }

            //Считаем определитель (произведение элементов на главной диагонали)
            det *= A[i][i];

            //Преобразуем последующие строки
            for (int k = i + 1; k < n; k++) {
                double factor = A[k][i] / A[i][i];
                //Для A нужно обновлять элементы от i до n, так как ранее элементы до i уже обнулились.
                for (int j = i; j < n; j++) {
                    A[k][j] -= A[i][j] * factor;
                }
                //Для identityMatrix обновляем всю строку (j от 0 до n)
                for (int j = 0; j < n; j++) {
                    identityMatrix[k][j] -= identityMatrix[i][j] * factor;
                }
                //Обновляем значения в столбце свободных членов
                b[k] -= b[i] * factor;
            }
        }

        //Сначала обратный ход метода Гаусса для поиска решения
        double[] answer = new double[n];
        answer[n - 1] = b[n - 1] / A[n - 1][n - 1]; //Решаем для последней переменной

        //Считаем последующие переменные
        for (int i = n - 2; i >= 0; i--) {
            answer[i] = b[i];
            for (int j = i + 1; j < n; j++) {
                answer[i] -= A[i][j] * answer[j];
            }
            answer[i] /= A[i][i];
        }

        System.arraycopy(answer, 0, solution, 0, n); //Копируем решение системы


        //Затем обратный ход метода Гаусса для поиска обратной матрицы (подставляем столбцы единичной матрицы вместо столбца свободных членов)
        for (int k = n - 1; k >= 0; k--) {
            answer[n - 1] = identityMatrix[n - 1][k] / A[n - 1][n - 1]; //Решаем для последней переменной

            //Считаем последующие переменные
            for (int i = n - 2; i >= 0; i--) {
                answer[i] = identityMatrix[i][k];
                for (int j = i + 1; j < n; j++) {
                    answer[i] -= A[i][j] * answer[j];
                }
                answer[i] /= A[i][i];
            }

            for (int i = 0; i < n; i++) inverseMatrix[i][k] = answer[i]; //Копируем весь столбец обратной матрицы
        }

        //Вычисление невязок
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += oldA[i][j] * solution[j];
            }
            residuals[i] = Math.abs(sum - oldB[i]); //Используем старые коэффициенты
        }

        //Учитываем количество перестановок строк для корректного знака определителя
        det = swapCount % 2 == 0 ? det : -det;
        return new Result(solution, inverseMatrix, det, residuals);
    }
}
