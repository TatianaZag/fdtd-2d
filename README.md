# fdtd-2d
В названии папок содержится наименование библиотек, которые были использованы для распараллеливания fdtd-2d. В корне каталога находится программа в ее первоначальном виде.

# Условные обозначения:
+ ✔️ — готово.
+ ⭕ — на стадии проверки.
+ ➕ — на стадии разработки.

# Реализация параллельной программы:
1. С помощью OpenMP:
    + ✔️ Классический метод с использованием директив parallel и for;
    + ✔️ С использованием механизма задач (директивы task, taskgroup и taskloop).
2. С помощью OpenACC:
    + ⭕ Версия для сервера с компилятором PGI и картой nvidia.

## OpenMP
Можно запустить с помощью Makefile в папке OpenMP (на данный момент запустится fdtd-2d-tasks.c, для запуска fdtd-2d-openmp.c следует сменить названия внутри Makefile на fdtd-2d-openmp.c):
```
make clean && make
```
Файл run.sh позволяет запустить fdtd-2d-tasks.c с разным числом задач на разных наборах данных и получить среднее время за десять запусков программы (можно поменять количество запусков в переменной ```N```):
```
chmod a+x run.sh
./run.sh
```
