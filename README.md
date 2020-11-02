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

Документация модуля с приведением используемых формул доступна 
[по этой ссылке](https://kononovarseniy.github.io/fompy/).

**Важно!** Основной системой единиц для расчётов является СГС (Гауссова система единиц).

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

В этом разделе представлены несколько примеров, демонстрирующих применение пакета FOMpy
для решения простейших типичных задач.

Подробное описание доступных подпрограмм, классов и их методов можно прочитать 
[по этой ссылке](https://kononovarseniy.github.io/fompy/).
В документации методов приведены используемые формулы и уравнения.

### Пример 1: расчёт концентрации дырок в легированном полупроводнике

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
  равную значению 1 эВ в единицах СГС.
  
- Находим концентрацию дырок (температура *T* = 300 К по умолчанию, уровень Ферми вычисляется автоматически):
  ```python
  >>> n_p = si_p.p_concentration()
  >>> print("{:e}".format(n_p))
  3.778950e+15
  ```

### Пример 2: определение проводимости

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

- Находим проводимость с помощью подпрограммы ```fompy.models.conductivity(n, n_mob, p, p_mob)```:
  ```python
  >>> sigma = conductivity(2. * 10**16, 3900. / volt, 8.5 * 10**16, 1900. / volt)
  >>> print("{:e}".format(sigma))
  3.448655e+13
  ```
  **Важно!** Так как в условии подвижность указана в единицах 
  см<sup>2</sup> В<sup>&minus;1</sup> с<sup>&minus;1</sup>, 
  то для перевода в СГС нужно домножить её на величину ```1. / fompy.constants.volt```, 
  равную значению 1 В<sup>&minus;1</sup> в системе СГС.

## Лицензия

[MIT](LICENSE.md)

## Помощь проекту и поддержка пользователей

Если вы желаете внести свой вклад в проект, следуйте инструкциям в файле [CONTRIBUTING.md](CONTRIBUTING.md). 
На репозитории действуют [правила поведения](CODE_OF_CONDUCT.md).

Предложения и пожелания по функционалу, наполнению проекта и исправлению ошибок принимаются на сайте репозитория 
в разделе [Issues](https://github.com/kononovarseniy/fompy/issues).
