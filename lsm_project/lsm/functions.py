"""
В этом модуле хранятся функции для применения МНК
"""

from math import sqrt

from typing import Optional
from numbers import Real       # раскомментируйте при необходимости

from lsm_project.event_logger.event_logger import EventLogger

from lsm_project.lsm.enumerations import MismatchStrategies
from lsm_project.lsm.models import (
    LSMDescription,
    LSMStatistics,
    LSMLines,
)


PRECISION = 3                   # константа для точности вывода
event_logger = EventLogger()    # для логирования


def get_lsm_description(
    abscissa: list[float], ordinates: list[float],
    mismatch_strategy: MismatchStrategies = MismatchStrategies.FALL
) -> LSMDescription:
    """
    Функции для получения описания рассчитаной зависимости

    :param: abscissa - значения абсцисс
    :param: ordinates - значение ординат
    :param: mismatch_strategy - стратегия обработки несовпадения

    :return: структура типа LSMDescription
    """

    global event_logger

    event_logger.info("Validating of data")

    try:
        abscissa = list(abscissa)
        ordinates = list(ordinates)
    except TypeError:
        event_logger.error("Type of the arguments is not a \
                           list (or convertible in list)")
        raise TypeError('An arguments must be a lists')

    if not _is_valid_measurments(abscissa) or \
       not _is_valid_measurments(ordinates):
        event_logger.error("There are nonreal values in input data")
        raise ValueError("Some values are not real")

    abscissa, ordinates = _process_mismatch(abscissa, ordinates,
                                            mismatch_strategy)

    event_logger.info("Validating passed successful")
    event_logger.info("Calculating coefficients started")

    data_lenght = len(abscissa)
    stats = _get_lsm_statistics(abscissa, ordinates)

    incline = (
        (stats.product_mean - stats.abscissa_mean * stats.ordinate_mean) /
        (stats.abs_squared_mean - stats.abscissa_mean**2)
    )

    shift = stats.ordinate_mean - incline * stats.abscissa_mean

    ordinates_error = (
        1 / (data_lenght - 2) *
        sum((ordinates[i] - incline * abscissa[i] - shift) ** 2
            for i in range(data_lenght))
    )

    incline_error = sqrt(ordinates_error / data_lenght /
                         (stats.abs_squared_mean - stats.abscissa_mean**2))
    shift_error = sqrt(ordinates_error * stats.abs_squared_mean / data_lenght /
                       (stats.abs_squared_mean - stats.abscissa_mean**2))

    event_logger.info("Calculating coefficients completed")

    return LSMDescription(
        incline=incline,
        shift=shift,
        incline_error=incline_error,
        shift_error=shift_error
    )


def get_lsm_lines(
    abscissa: list[float], ordinates: list[float],
    lsm_description: Optional[LSMDescription] = None
) -> LSMLines:
    """
    Функция для расчета значений функций с помощью результатов МНК

    :param: abscissa - значения абсцисс
    :param: ordinates - значение ординат
    :param: lsm_description - описание МНК

    :return: структура типа LSMLines
    """
    if lsm_description is None:
        lsm_description = get_lsm_description(abscissa, ordinates)

    if not isinstance(lsm_description, LSMDescription):
        event_logger.error('lms_descriptor takes object of unannotated type')
        raise TypeError('Incorrect data in lms_description')

    line_predicted = [lsm_description.incline * x + lsm_description.shift
                      for x in abscissa]
    line_above = [(lsm_description.incline + lsm_description.incline_error) *
                  x + lsm_description.shift + lsm_description.shift_error
                  for x in abscissa]
    line_under = [(lsm_description.incline - lsm_description.incline_error) *
                  x + lsm_description.shift - lsm_description.shift_error
                  for x in abscissa]

    return LSMLines(
        abscissa=abscissa,
        ordinates=ordinates,
        line_predicted=line_predicted,
        line_above=line_above,
        line_under=line_under
    )


def get_report(
    lsm_description: LSMDescription, path_to_save: str = None
) -> str:
    """
    Функция для формирования отчета о результатах МНК

    :param: lsm_description - описание МНК
    :param: path_to_save - путь к файлу для сохранения отчета

    :return: строка - отчет определенного формата
    """
    global PRECISION

    report = '\n'.join([
        'LSM computing result'.center(100, "="), '',
        f"[INFO]: incline: {lsm_description.incline:.{PRECISION}f};",
        f"[INFO]: shift: {lsm_description.shift:.{PRECISION}f};",
        f"[INFO]: incline error: {lsm_description.incline_error:.{PRECISION}f};",
        f"[INFO]: shift error: {lsm_description.shift_error:.{PRECISION}f};",
        '', ''.center(100, "=")
    ])

    if path_to_save:
        with open(path_to_save, 'w') as fts:
            fts.write(report)

    return report


# служебная функция для валидации
def _is_valid_measurments(measurments: list[float]) -> bool:
    return all(isinstance(el, Real) for el in measurments)


# служебная функция для обработки несоответствия размеров
def _process_mismatch(
    abscissa: list[float], ordinates: list[float],
    mismatch_strategy: MismatchStrategies = MismatchStrategies.FALL
) -> tuple[list[float], list[float]]:
    global event_logger

    if len(abscissa) < 3 or len(ordinates) < 3:
        event_logger.error("A lenght of arguments is less than 3")
        raise ValueError("Arguments lenghts must be more than 2")

    if len(abscissa) != len(ordinates):
        match mismatch_strategy:
            case MismatchStrategies.FALL:
                event_logger.error("Lenghts of the arguments are not equal")
                raise RuntimeError("Lenghts of the arguments are not equal")
            case MismatchStrategies.CUT:
                event_logger.warning("Lenghts of the arguments are not equal, \
                                     so longest list will be cut")
                min_lenght = min(map(len, (abscissa, ordinates)))
                abscissa, ordinates = (abscissa[:min_lenght],
                                       ordinates[:min_lenght])
            case _:
                event_logger.error('Unknown data in argument "mismatch_data"')
                raise ValueError('Unknown data in argument "mismatch_data"')

    return abscissa, ordinates


# служебная функция для получения статистик
def _get_lsm_statistics(
    abscissa: list[float], ordinates: list[float]
) -> LSMStatistics:
    global event_logger, PRECISION

    data_lenght = len(abscissa)

    abscissa_mean = sum(abscissa) / data_lenght
    ordinates_mean = sum(ordinates) / data_lenght
    product_mean = sum(abscissa[i] * ordinates[i]
                       for i in range(data_lenght)) / data_lenght
    abs_squared_mean = sum(abscissa[i]**2
                           for i in range(data_lenght)) / data_lenght

    return LSMStatistics(
        abscissa_mean=abscissa_mean,
        ordinate_mean=ordinates_mean,
        product_mean=product_mean,
        abs_squared_mean=abs_squared_mean
    )


# служебная функция для получения описания МНК
def _get_lsm_description(
    abscissa: list[float], ordinates: list[float]
) -> LSMDescription:
    global event_logger, PRECISION

    # ваш код
    # эту строчку можно менять
    return LSMDescription(
        incline=0,
        shift=0,
        incline_error=0,
        shift_error=0
    )
