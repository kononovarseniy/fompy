[![PyPI version](https://img.shields.io/pypi/v/fti-fompy.svg)](https://pypi.python.org/pypi/fti-fompy/)
[![PyPI license](https://img.shields.io/pypi/l/fti-fompy.svg)](https://pypi.python.org/pypi/fti-fompy/)
![GitHub contributors](https://img.shields.io/github/contributors/kononovarseniy/fompy)<br/>
![PyPI - Downloads](https://img.shields.io/pypi/dm/fti-fompy?label=pypi%20downloads)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/kononovarseniy/fompy)

# FOMpy
FOMpy - подпрограммы для курса "Физические Основы Микроэлектроники"
Предполагается что мы совместными усилиями сможем создать достаточную базу подпрограмм,
для решения задач по ФОМЭ. Важным отличием данного проекта является использование СГС как
основной системы единиц.

Описание доступных классов и их методов можно прочитать [здесь](https://kononovarseniy.github.io/fompy/).

Если вы желаете внести свой вклад в проект, следуйте инструкциям в файле [CONTRIBUTING.md](CONTRIBUTING.md).

## Установка
### Глобальная установка
Следуйте этой инструкции если хотите использовать FOMpy в своём проекте или как калькулятор.

Если вы хотите просто установить FOMpy глобально достаточно выполнить следующую команду:
```
$ pip install fti-fompy
```
Однако, рекомендуется использовать виртуальную среду, чтобы избежать возможных проблем с установкой зависимостей.  

### Установка в виртуальной среде
- Создайте виртуальную среду
    ```
    $ python -m venv .venv
    ```
    Возможно вам необходимо будет написать ```python3``` вместо ```python```
- Запустите виртуальную среду:
    ```
    $ source ./venv/bin/activate
    ```
    Эта команда активирует виртуальную среду. Её нужно выполнять каждый раз перед запуском скриптов.
    Эффект действует до закрытия окна терминала, или вызова команды ```deactivate```
- Установите FOMpy:
    ```
    $ pip install fti-fompy
    ```

### Удобный скрипт для запуска
- Рекомендую добавить в свой файл ```~/.bashrc``` следующие строки
    ```
    FOMPY_IMPORTS="
    from math import *
    from fompy.constants import *
    from fompy.materials import *
    from fompy.models import *
    from fompy.units import unit
    "
    
    fompy() {
        cd <Путь до папки с FOMpy> # Эти две строки нужны только для 
        source .venv/bin/activate  # запуска виртуальной среды
        PYTHONSTARTUP=<(echo "$FOMPY_IMPORTS") python
    }
    ```
    Это позволит использовать в терминале команду ```fompy``` для запуска настроенного интерпретатора python.
    Для применения изменений перезапустите терминал.

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
