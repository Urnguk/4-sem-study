# Молекулярное моделирование Аргона
Цель проекта: проверить выполнение термодинамических соотношений, используя метод молекулярного моделирования и законы Ньютоновской механики.

Команда разработчиков: Горшков Тимофей, Потапов Денис, Кутушева Алиса, Лапа Ксения.

Список использованных статей и литературы: 
  Michiel Bosch (2007). The Molecular Dynamic Simulation of neutral Argon Particles Michiel Bosch (2007) 
  Verlet, Loup (1967). "Computer "Experiments" on Classical Fluids 
  Victor Ruhle (August 8, 2007). Berendsen and Nose-Hoover thermostats 
  Д.В. Сивухин. Общий курс физики. Термодинамика и молекулярная физика. т. II. 
  https://en.wikipedia.org/wiki/Periodic_boundary_conditions

Инструкция по запуску:
  Открыть файл Argon.cpp в MS Visual Studio 2015 и скомпилировать в Release x86.

План выполнения: 
1) Создать рабочую модель реального газа. // Done
2) Добиться выполнения ЗСЭ в системе. // Done
3) Привести систему в квазиравновесное состояние (используя термостат), для проверки термодинамических соотношений. // In process (Потапов Денис) 10.05.18
4) Проверка соответствия распределения по скоростям в системе распределению Максвелла. // Done
5) Проверка соответствия случайных блужданий молекулы закону диффузии Эйнштейна-Смолуховского. // Done
6) Визуализация полученных данных в виде графиков. // In process (Горшков Тимофей) 10.05.18
7) Оптимизация рабочего кода. // In process (Горшков Тимофей, Потапов Денис) 17.05.18
7) Визуализация траекторий движения молекул в объеме. // In perspective 
8) Создание графического интерфейса. // In perspective
