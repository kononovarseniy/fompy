# FOMpy
FOMpy - подпрограммы для курса "Физические Основы Микроэлектроники"
Предполагается что мы совместными усилиями сможем создать достаточную базу подпрограмм,
для решения задач по ФОМЭ.

На данный момент предпологается использование пакета через терминал.
Но любители Jupiter notebook могут предоставить инструкции для установки и использования. 

## Правила репозитория
- Все величины и расчёты в СГС (Гауссовой системе единиц),
    даже в промежуточных вычисления применение внесистемных единиц не приветствуется.
- Все величины физических констант приводить с точностью до 4 значащих цифр.
- Названия функций на английском.
- Документация к каждой функции (какую величину вычисляет, какие принимает),
    можно на русском.
- Желательно подписывать каждый коммит используя опцию -S

## Использование через терминал
### Установка
После клонирования репозитория, для использования скриптов необходимо,
выполнить слудующие команды в директории проекта.
- Создать виртуальную среду для установки пакетов:
    ```
    python -m venv .venv
    ```
  Возможно вам необходимо будет написать ```python3``` вместо ```python```
- Запустить виртуальную среду:
    ```
    source ./venv/bin/activate
    ```
- Уствновить зависимости:
    ```
    pip install --upgrade pip wheel
    pip install -r requirements.txt
    ```
### Запуск
- Если вы закрыли терминал или деактивировали виртуальную среду (комманда ```deactivate```),
    повторите второй шаг процедуры установки:
    ```
    source ./venv/bin/activate
    ```
- Теперь можно запустить интерпретатор питона, чтобы использовать его как калькулятор:
    ```
    python
    >>> from math import * # Определяет свою константу 'e' поэтому импортируется раньше остальных
    >>> from constants import *
    >>> from phys import *
    >>> E = me * c**2
    ```

### Расчет параметров полупроводников
В модуле ```materials``` есть определения нескольких материалов (пока только Si).
Новый материал можно определеить создав экземпляр класса ```materials.Semiconductor```
```
Semiconductor(electron_effective_mass. hole_effective_mass, energy_gap, electron_affinity)
```
Если параметры вам неизвестны (и не участвуют в расчётах) вы можете передать None в качестве этого параметра.

Класс полупроводника объявляет несколько методов:
- ```Nc(self, T=300)``` -- Еффективная плотность состояний в зоне проводимости.
- ```Nv(self, T=300)``` -- Еффективная плотность состояний в валентной зоне.
- ```p_intrinsic(self, Ef=None, T=300)``` -- собсвенная концентрация дырок.
- ```n_intrinsic(self, Ef=None, T=300)``` -- собственная концентрауий електронов.
- ```fermi_level(self, T=300)``` -- уровень ферми из условия электронейтральности.

Большинство этих методов принимают уровень ферми ```Ef``` и температуру ```T```.
Оба аргумента необязательны, по умолчанию ```T=300```, а ```Ef``` вычисляется методом ```fermi_level```.

Класс примесного полупроводника ```materials.DopedSemiconductor``` расширяет класс проводника и добавляет следующие методы:
- ```DopedSemiconductor(self, mat, Na, Ea, Nd, Ed)``` -- Конструктор класса.
- ```p_donor_concentration(self, Ef=None, T=300)``` -- Концентрация положительных ионов доноров.
- ```n_acceptor_concentration(self, Ef=None, T=300)``` -- Концентрация отрицательных ионов аццепторов.
