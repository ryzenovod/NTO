import sys

def solve():
    """
    Решение задачи "Непрерывные подтрезки" с использованием алгоритма скользящего окна.
    """
    try:
        # Быстрое чтение ввода для больших данных
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str)
        a = list(map(int, sys.stdin.readline().split()))
    except (IOError, ValueError):
        # Обработка пустого ввода или ошибок
        return

    # Шаг 1: Вычисляем префиксные суммы для быстрого нахождения суммы любого подмассива.
    # prefix_sums[i] хранит сумму элементов a[0]...a[i-1].
    prefix_sums = [0] * (n + 1)
    for i in range(n):
        prefix_sums[i + 1] = prefix_sums[i] + a[i]

    # Шаг 2: Инициализируем переменные для алгоритма скользящего окна.
    max_len = 0          # Максимальная длина найденного уникального подмассива.
    max_sum = 0          # Максимальная сумма для подмассива длины max_len.
    left = 0             # Левый указатель окна.
    last_seen = {}       # Словарь для хранения последних индексов чисел: {число: индекс}.

    # Шаг 3: Проходим по массиву правым указателем, расширяя окно.
    for right in range(n):
        current_element = a[right]

        # Если текущий элемент уже встречался в пределах нашего окна (т.е. его последний
        # индекс >= left), то мы нашли дубликат.
        if current_element in last_seen and last_seen[current_element] >= left:
            # Сдвигаем левую границу окна вправо за предыдущее вхождение этого элемента.
            left = last_seen[current_element] + 1
        
        # Обновляем индекс последнего вхождения текущего элемента.
        last_seen[current_element] = right
        
        # Теперь окно a[left...right] гарантированно содержит уникальные элементы.
        # Вычисляем его длину и сумму.
        current_len = right - left + 1
        current_sum = prefix_sums[right + 1] - prefix_sums[left]

        # Шаг 4: Обновляем результат в соответствии с условиями задачи.
        # Приоритет — длина. Сумма — вторичный критерий.
        if current_len > max_len:
            max_len = current_len
            max_sum = current_sum
        elif current_len == max_len:
            # Если нашли подмассив такой же максимальной длины, выбираем тот, у которого сумма больше.
            max_sum = max(max_sum, current_sum)

    # Шаг 5: Согласно условию, подмассив должен состоять не менее чем из двух элементов.
    if max_len < 2:
        print(0)
    else:
        print(max_sum)

solve()