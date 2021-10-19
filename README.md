[![PyPI version](https://img.shields.io/pypi/v/fti-fompy.svg)](https://pypi.python.org/pypi/fti-fompy/)
[![PyPI license](https://img.shields.io/pypi/l/fti-fompy.svg)](https://pypi.python.org/pypi/fti-fompy/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/fti-fompy?label=pypi%20downloads) <br/>
[![GitHub contributors](https://img.shields.io/github/contributors/kononovarseniy/fompy)](https://github.com/kononovarseniy/fompy/graphs/contributors)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/kononovarseniy/fompy)
![Coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fgist.githubusercontent.com%2Fkononovarseniy%2F44e4ca5d46404d5c37ab1b8661bd6675%2Fraw%2Fcoverage.json) <br/>
[![GitHub milestone](https://img.shields.io/github/milestones/progress/kononovarseniy/fompy/1)](https://github.com/kononovarseniy/fompy/milestone/1)

# FOMpy

FOMpy — подпрограммы и классы для курса «Физические основы микроэлектроники» (ФОМЭ).
Идея проекта в том, чтобы совместными усилиями создать достаточную базу подпрограмм для решения задач по ФОМЭ.

Документация модуля с приведением используемых формул и единиц измерения доступна 
[по этой ссылке](https://kononovarseniy.github.io/fompy/).

**Важно!** Основной системой единиц для расчётов является СГС (Гауссова система единиц).
О том, как работать с единицами измерения, читайте в разделе [Единицы измерения](#единицы-измерения).

## Установка

Данная инструкция поможет вам начать использовать пакет FOMpy в вашем проекте 
или работать с ним как с калькулятором.

С целью избежать возможных проблем с установкой зависимостей,
рекомендуется [установка в виртуальной среде](#установка-в-виртуальной-среде).

### Глобальная установка

Для того чтобы установить FOMpy глобально, выполните команду
```console
$ pip install fti-fompy
```

### Установка в виртуальной среде

- Создайте виртуальную среду (возможно, придётся написать ```python3``` вместо ```python```):
    ```console
    $ python -m venv .venv
    ```

- Запустите виртуальную среду:
    ```console
    $ source ./venv/bin/activate
    ```
  **Важно!** В дальнейшем эту команду нужно будет выполнять каждый раз перед запуском скриптов, работающих с FOMpy.
  Эффект действует до закрытия окна терминала или вызова команды ```deactivate```.
  
- Установите пакет FOMpy:
    ```console
    $ pip install fti-fompy
    ```

### Удобный скрипт для запуска

Вы можете настроить терминал таким образом, чтобы одной командой в нём 
запускался интерпретатор Python, сразу после запуска готовый к работе с FOMpy.

- Добавьте в файл ```~/.bashrc``` (или другой rc-файл в зависимости от вашей командной оболочки) следующие строки:
    ```sh
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
  
- Перезапустите терминал.

- Наберите команду
    ```console
    $ fompy
    ```
  Теперь этой командой вы можете вызывать интерпретатор ```python```, 
  в котором уже будут импортированы все нужные модули FOMpy.

## Использование

В этом разделе объясняется, как наиболее эффективно работать с единицами измерения, 
а также представлены несколько примеров, демонстрирующих применение пакета FOMpy
для решения простейших типичных задач.

Подробное описание доступных подпрограмм, классов и их методов можно прочитать 
[по этой ссылке](https://kononovarseniy.github.io/fompy/).
В документации методов приведены используемые формулы и уравнения, 
а также указаны единицы измерения физических величин.

### Единицы измерения

С единицами физических величин в расчётах необходимо обращаться осторожно. 
Важно помнить, что в пакете FOMpy **<u>везде используется система СГС</u>**: 
как в возвращаемых значениях функций, так и в аргументах, передаваемых тем же функциям. 
Все физические константы также приведены в единицах СГС.

Функционал, облегчающий работу с единицами измерения, расположен 
в модулях `fompy.constants` и `fompy.units`. 
Документация методов и подпрограмм снабжена обозначениями принимаемых и возвращаемых единиц 
в квадратных скобках, например, 
\[g\] &mdash; грамм, \[statV\] &mdash; статвольт, \[1\] &mdash; безразмерная единица.

#### Функция `fompy.units.unit()`

Общий интерфейс для работы с единицами предоставляется функцией `fompy.units.unit()`. 
Пример использования:
```python
s = Semiconductor(me_eff= 0.1 * unit('MeV_m'), mh_eff= 4. * 10**(-28), Eg= 1. * unit('eV'))
```
Здесь масса электрона задана эквивалентным значением энергии покоя в МэВ, 
масса дырки сразу приведена в граммах, а размер запрещённой зоны равен 1 эВ.

Строка, передаваемая функции `unit()`, может содержать
* наименование единицы измерения (обязательно);
* дольную или кратную приставку в виде префикса к наименованию;
* рациональная степень;
* пробельные символы (где угодно, но **не** между приставкой и наименованием единицы).

Показатель степени по умолчанию равен 1, в иных случаях он присоединяется в виде суффикса справа от 
наименования единицы. Этот суффикс выглядит так: `^` + знак (`-` для отрицательной степени, 
`+` для положительной степени) + целое или рациональное число. Знаки умножения, возведения в степень, и знак `+` в показателе степени &mdash; необязательны.

Примеры: микрометр в минус третьей степени &mdash; `um^-3` или `um-3`; вольт &mdash; `V`, или `V ^1`, или `V ^ +1`; вольт в степени 3/2 &mdash; `V ^ 3/2`, `V 3/2`.

Комбинации разных единиц строятся с помощью символов `/` для частного и `*`/отсутствие символа для произведения.
Важный момент: символ `/` делит обозначение на числитель и знаменатель. Всё, что записано после него, станет множителем с 
противоположной по знаку степенью. В качестве числителя может выступать `1`.

Функция всегда способна отличить знак `/` в показателе степени от знака `/` разделяющего числитель и знаменатель, поскольку знаменатель не может начинаться с числа.

Примеры: ом-метр &mdash; `Ohm m`, или `Ohm * m`, или `Ohm / m-1`; 
секунда в минус первой степени &mdash; `s^-1` или `1 / s`; 
генри = кг м<sup>2</sup> с<sup>&minus;2</sup> А<sup>&minus;2</sup> &mdash; `kg m^2 / s^2 A^2`.

Поддерживаются следующие наименования единиц измерения:

Подстрока | Единица
--- | --- 
`m` | метр
`g` | грамм
`s` | секунда
`K` | кельвин
`A` | ампер
`eV` | электронвольт
`eV_m` | 1 эВ/c<sup>2</sup> (масса, эквивалентная энергии покоя 1 эВ)
`eV_T` | 1 эВ/k<sub>B</sub> (температура, эквивалентная тепловой энергии 1 эВ)
`Hz` | герц
`N` | ньютон
`J` | джоуль
`W` | ватт
`Pa` | паскаль
`C` | кулон
`V` | вольт
`Ohm` | ом
`F` | фарад
`Wb` | вебер
`T` | тесла
`H` | генри

Таблица дольных и кратных десятичных приставок:

Подстрока | Множитель | Приставка
--- | --- | ---
`y` | 10<sup>&minus;24</sup> | иокто-
`z` | 10<sup>&minus;21</sup> | зепто-
`a` | 10<sup>&minus;18</sup> | атто-
`f` | 10<sup>&minus;15</sup> | фемто-
`p` | 10<sup>&minus;12</sup> | пико-
`n` | 10<sup>&minus;9</sup> | нано-
`u` | 10<sup>&minus;6</sup> | микро-
`m` | 10<sup>&minus;3</sup> | милли-
`c` | 10<sup>&minus;2</sup> | санти-
`d` | 10<sup>&minus;1</sup> | деци-
`da` | 10<sup>1</sup> | дека-
`h` | 10<sup>2</sup> | гекто-
`k` | 10<sup>3</sup> | кило-
`M` | 10<sup>6</sup> | мега-
`G` | 10<sup>9</sup> | гига-
`T` | 10<sup>12</sup> | тера-
`P` | 10<sup>15</sup> | пета-
`E` | 10<sup>18</sup> | экса-
`Z` | 10<sup>21</sup> | зетта-
`Y` | 10<sup>24</sup> | иотта-

#### Физические константы

Физические константы определены в модуле `fompy.constants`. Все значения приведены в единицах системы СГС.

Таблица представленных констант:

Имя объекта | Константа
--- | ---
`c` | скорость света
`k` | постоянная Больцмана
`Na` | постоянная Авогадро
`R` | универсальная газовая постоянная
`sigma` | константа Стефана-Больцмана
`me` | масса электрона
`mp` | масса протона
`mn` | масса нейтрона
`e` | заряд электрона
`h` | постоянная Планка
`h_bar` | редуцированная постоянная Планка
`eV` | 1 эВ
`eV_m` | 1 эВ/c<sup>2</sup> (масса, эквивалентная энергии покоя 1 эВ)
`eV_T` | 1 эВ/k<sub>B</sub> (температура, эквивалентная тепловой энергии 1 эВ)
`amu` | атомная единица массы (дальтон)
`angstrom` | ангстрем
`volt` | 1 В
`ampere` | 1 А
`ohm` | 1 Ом
`farad` | 1 Ф
`henry` | 1 Гн
`Ry` | постоянная Ридберга
`a0` | радиус Бора

### Примеры

#### Пример 1: расчёт концентрации дырок в легированном полупроводнике

Необходимо найти концентрацию дырок в кремнии Si, легированном акцепторной примесью.
Температура *T* = 300 К, концентрация акцепторной примеси 
*N<sub>a</sub>* = 10<sup>17</sup> см<sup>&minus;3</sup>,
акцепторный уровень *E<sub>a</sub>* = 0,3 эВ (от вершины валентной зоны *E<sub>v</sub>* = 0).

- Подготавливаем интерпретатор Python — воспользуемся [удобной настройкой Unix shell](#удобный-скрипт-для-запуска):
  ```console
  $ fompy
  ```

- Создаём объект легированного полупроводника ```fompy.models.DopedSemiconductor```: 
  * задаём в качестве базового материала кремний ```fompy.materials.Si```;
  * приводим значения для акцепторных концентрации *N<sub>a</sub>* и уровня *E<sub>a</sub>*; 
  * зануляем параметры донорной примеси *N<sub>d</sub>* и *E<sub>d</sub>* (считаем, что она отсутствует).
  ```python
  >>> si_p = DopedSemiconductor(Si, 10**17, 0.3 * eV, 0, 0)
  ```
  **Важно!** Применение эВ требует домножения на величину ```fompy.constants.eV```, 
  равную значению 1 эВ в единицах СГС. Либо можно взамен воспользоваться функцией ```fompy.units.unit()```, описанной 
  в разделе [Единицы измерения](#единицы-измерения).
  
- Находим концентрацию дырок (температура *T* = 300 К по умолчанию, уровень Ферми вычисляется автоматически):
  ```python
  >>> n_p = si_p.p_concentration()
  >>> print("{:e}".format(n_p))
  3.778950e+15
  ```

#### Пример 2: определение проводимости

Требуется вычислить проводимость материала при заданных концентрации 
электронов *n<sub>n</sub>* = 2,0 &middot; 10<sup>16</sup> см<sup>&minus;3</sup> 
и дырок *n<sub>p</sub>* = 8,5 &middot; 10<sup>16</sup> см<sup>&minus;3</sup>, а также
подвижности электронов &mu;*<sub>n</sub>* = 3,9 &middot; 10<sup>3</sup> 
см<sup>2</sup> В<sup>&minus;1</sup> с<sup>&minus;1</sup> 
и дырок &mu;*<sub>p</sub>* = 1,9 &middot; 10<sup>3</sup> 
см<sup>2</sup> В<sup>&minus;1</sup> с<sup>&minus;1</sup>.

- Подготавливаем интерпретатор Python:
  ```console
  $ fompy
  ```

- Находим проводимость с помощью подпрограммы ```fompy.models.conductivity(n, n_mob, p, p_mob)```. Здесь для перевода 
  [единиц измерения](#единицы-измерения) мы пользуемся функцией ```fompy.units.unit()```.
  ```python
  >>> sigma = conductivity(2. * 10**16, 3900. * unit('cm^2 / V s'), 8.5 * 10**16, 1900. * unit('cm^2 / V s'))
  >>> print("{:e}".format(sigma))
  3.448655e+13
  ```

## Лицензия

[MIT](LICENSE.md)

## Помощь проекту и поддержка пользователей

Если вы желаете внести свой вклад в проект, следуйте инструкциям в файле [CONTRIBUTING.md](CONTRIBUTING.md). 
На репозитории действуют [правила поведения](CODE_OF_CONDUCT.md).

Предложения и пожелания по функционалу, наполнению проекта и исправлению ошибок принимаются на сайте репозитория 
в разделе [Issues](https://github.com/kononovarseniy/fompy/issues).
