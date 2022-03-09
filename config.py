import re
import pyfiglet


def clear_config():
    res = pyfiglet.figlet_format('approximate microscopic model of two-dim. filtration\n Config file', font="slant",
                                 justify='center')
    with open('config.dat', 'w') as file:
        file.write(res)
        file.write("""
N = 10
M = 10
gu_const_val('right','p', 0)
gu_const_val('top','q', 0)
gu_const_val('left','p', 10)
gu_const_val('bottom','q', 0)""")
    print('Конфигурационный файл сброшен')


def read_config():
    global N, M
    global right_func, right_values
    global top_func, top_values
    global left_func, left_values
    global bottom_func, bottom_values
    try:
        with open('config.dat', 'r') as file:
            for i, line in enumerate(file):
                if i == 37:
                    N = int(re.findall(r'\d+', line)[0])
                elif i == 38:
                    M = int(re.findall(r'\d+', line)[0])
                elif i == 39:
                    right_func = re.findall(r'gu_\w+', line)[0]
                    right_values = re.sub(r"'", r'', re.findall(r'\(.+\)', line)[0][1:-1]).split(',')
                elif i == 40:
                    top_func = re.findall(r'gu_\w+', line)[0]
                    top_values = re.sub(r"'", r'', re.findall(r'\(.+\)', line)[0][1:-1]).split(',')
                elif i == 41:
                    left_func = re.findall(r'gu_\w+', line)[0]
                    left_values = re.sub(r"'", r'', re.findall(r'\(.+\)', line)[0][1:-1]).split(',')
                elif i == 42:
                    bottom_func = re.findall(r'gu_\w+', line)[0]
                    bottom_values = re.sub(r"'", r'', re.findall(r'\(.+\)', line)[0][1:-1]).split(',')
            else:
                if i < 42:
                    raise Exception('Некорректный конфигурационный файл')
    except Exception as ex:
        print('Ошибка при чтении конфигурационного файла:', ex)
        clear_config()
        read_config()


if __name__ == "__main__":
    read_config()
    print(N, M)
    print(right_func, right_values)
    print(top_func, top_values)
    print(left_func, left_values)
    print(bottom_func, bottom_values)
