[![PyPI version shields.io](https://img.shields.io/pypi/v/fti-fompy.svg)](https://pypi.python.org/pypi/fti-fompy/)
[![PyPI license](https://img.shields.io/pypi/l/fti-fompy.svg)](https://pypi.python.org/pypi/fti-fompy/)

# FOMpy
FOMpy - подпрограммы для курса "Физические Основы Микроэлектроники"
Предполагается что мы совместными усилиями сможем создать достаточную базу подпрограмм,
для решения задач по ФОМЭ.

## Правила репозитория
- Все величины и расчёты в СГС (Гауссовой системе единиц),
    даже в промежуточных вычисления применение внесистемных единиц не приветствуется.
- Все величины физических констант приводить с точностью до 4 значащих цифр.
- Названия функций на английском.
- Документация к каждой функции (какую величину вычисляет, какие принимает),
    можно на русском.
- Желательно подписывать каждый коммит используя опцию -S

## Установка (Для использования FOMpy)
Следуйте этой инструкции если хотите использовать FOMpy в своём проекте или как калькулятор.
- Создать виртуальную среду (можете пропустить этот шаг если вы хотите установить пакет глобально)
    ```
    $ python -m venv .venv
    ```
    Возможно вам необходимо будет написать ```python3``` вместо ```python```
- Запустить виртуальную среду:
    ```
    $ source ./venv/bin/activate
    ```
    Эта команда активирует виртуальную среду. Её нужно выполнять каждый раз перед запуском скриптов.
    Эффект действует до закрытия окна терминала, или вызова команды ```deactivate```
- Установить FOMpy:
    ```
    $ pip install fti-fompy
    ```
### Удобный скрипт для запуска
- Рекомендую добавить в свой файл ```~/.bashrc``` следующие строки.
    Это позволит использовать в терминале команду ```fompy``` для запуска настроенного интерпретатора python.
    Для применения изменений перезапустите терминал.
    ```
    fompy() {
        cd <Полный путь до папки содержащей FOMpy>
        source .venv/bin/activate
        PYTHONSTARTUP=<(echo -e 'from math import *\nfrom fompy.materials import *\nfrom fompy.constants import *\nfrom fompy.phys import *\nfrom fompy.units import unit') python
    }
    ```

## Установка (Для разработки и тестирования FOMpy)
Это инструкция для тех кто хочет внести свой вклад в разработку FOMpy.
- Сначала клонируйте репозиторий используя команду
    ```
    $ git clone https://github.com/kononovarseniy/fompy.git
    ```
- Перейдите в директорию проекта:
    ```
    $ cd fompy
    ```
- Теперь надо создать виртуальную среду чтобы устанавливать зависимости локально:
    ```
    $ python -m venv .venv
    $ source ./venv/bin/activate
    $ pip install --upgrade pip wheel
    ```
- И установить все зависимости
    ```
    $ pip install -r requirements.txt
    ```
- Теперь можно открыть проект в вашей любимой IDE, например, в PyCharm.

## Применение FOMpy для расчета параметров полупроводников
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
